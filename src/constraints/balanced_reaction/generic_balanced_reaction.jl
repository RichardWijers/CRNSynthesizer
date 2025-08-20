struct LocalGenericBalancedReaction <: AbstractLocalConstraint
    path::Vector{Int}  # Path to the reaction node
    rule_to_dict::Dict{Int, Dict{String, Int}}  # Mapping of rule indices to atom counts
    hole_terminals::Vector{Int}  # Indices of terminal holes in the grammar
end

import Base.==
function ==(a::LocalGenericBalancedReaction, b::LocalGenericBalancedReaction)
    return a.path == b.path
end

import Base.hash
function hash(a::LocalGenericBalancedReaction, h::UInt)
    return hash(a.path, h)
end

function HerbConstraints.shouldschedule(
        solver::Solver, constraint::LocalGenericBalancedReaction, path::Vector{Int}
)
    # Check if the path is a child of the reaction node
    if length(path) <= length(constraint.path)
        return false
    end

    if path[1:length(constraint.path)] != constraint.path
        return false
    end

    # Schedule for atom nodes and relevant holes
    node = get_node_at_location(solver, path)
    type = get_node_type(solver.grammar, node)

    # Only schedule for holes that could contain atoms
    return type == :atom ||
           type == :chain ||
           type == :structure ||
           type == :branch ||
           type == :branches ||
           type == :molecule_list ||
           type == :molecule
end

mutable struct RelevantCounts
    atom_holes::Int
    atoms_fixed::Dict{String, Int}
    holes::Int
    hole_paths::Vector{Vector{Int}}
end

import Base: +
function +(a::RelevantCounts, b::RelevantCounts)::RelevantCounts
    # Merge atom count dicts by summing values for matching keys
    atoms = Dict{String, Int}()
    for (k, v) in a.atoms_fixed
        atoms[k] = get(atoms, k, 0) + v
    end
    for (k, v) in b.atoms_fixed
        atoms[k] = get(atoms, k, 0) + v
    end

    return RelevantCounts(
        a.atom_holes + b.atom_holes,
        atoms,
        a.holes + b.holes,
        vcat(a.hole_paths, b.hole_paths)
    )
end

function get_relevant_counts(
        solver::Solver, c::LocalGenericBalancedReaction, path::Vector{Int}, node::AbstractRuleNode, counts::RelevantCounts=RelevantCounts(0, Dict{String, Int}(), 0, Vector{Vector{Int}}())
)::RelevantCounts
    # Get the relevant paths for the current node
    node_type = get_node_type(solver.grammar, node)

    get_elem_from_atom_rule(rule_val) = begin
        s = String(rule_val)
        if startswith(s, "[") && endswith(s, "]")
            return s[2:(end - 1)]
        else
            return s
        end
    end

    if isuniform(node)
        @match node_type begin
            :molecule_list => begin
                rule = HerbConstraints.get_rule(node)
                @match solver.grammar.rules[rule] begin
                    :(vcat([molecule], molecule_list)) ||
                    :(vcat([required_molecule], molecule_list)) => begin
                        get_relevant_counts(
                            solver, c, child(path, 1), get_children(node)[1], counts
                        )
                        get_relevant_counts(
                            solver, c, push!(path, 2), get_children(node)[2], counts
                        )

                        return counts
                    end
                    :(Vector{Molecule}()) => begin
                        return counts
                    end
                    x => throw("Unknown molecule_list rule: $x")
                end
            end
            :molecule => begin
                get_relevant_counts(solver, c, push!(path, 1), node.children[1], counts)

                return counts
            end
            :chain => begin
                rule = HerbCore.get_rule(node)
                @match solver.grammar.rules[rule] begin
                    :(structure * bond * chain) => begin
                        # Get the relevant paths for the structure and chain
                        get_relevant_counts(
                            solver, c, child(path, 1), node.children[1], counts
                        )
                        get_relevant_counts(
                            solver, c, push!(path, 3), node.children[3], counts
                        )

                        return counts
                    end
                    :(SMILES_combine_chain(bond, structure, chain)) => begin
                        # Get the relevant paths for the bond and structure
                        get_relevant_counts(
                            solver, c, child(path, 2), node.children[2], counts
                        )
                        get_relevant_counts(
                            solver, c, push!(path, 3), node.children[3], counts
                        )

                        return counts
                    end

                    :(atom * ringbonds) => begin
                        # Get the relevant path for the atom
                        atom_path = push!(path, 1)
                        atom = get_node_at_location(solver, atom_path)

                        if atom isa RuleNode
                            # If the atom is filled, count it as a fixed atom
                            rule_idx = HerbCore.get_rule(atom)
                            rule_val = solver.grammar.rules[rule_idx]
                            elem = get_elem_from_atom_rule(rule_val)

                            counts.atoms_fixed[elem] = get(
                                counts.atoms_fixed, elem, 0
                            ) + 1

                            return counts
                        else
                            counts.atom_holes += 1
                            return counts
                        end
                    end
                    x => throw("Unknown chain rule: $x")
                end
            end
            :structure => begin
                get_relevant_counts(
                    solver, c, child(path, 3), node.children[3], counts
                )

                atom_path = push!(path, 1)
                atom = get_node_at_location(solver, atom_path)

                if isfilled(atom)
                    # If the atom is filled, count it as a fixed atom
                    rule_idx = HerbCore.get_rule(atom)
                    rule_val = solver.grammar.rules[rule_idx]
                    elem = get_elem_from_atom_rule(rule_val)

                    counts.atoms_fixed[elem] = get(
                        counts.atoms_fixed, elem, 0
                    ) + 1
                else
                    counts.atom_holes += 1
                end

                return counts
            end
            :ringbonds => begin
                return counts
            end
            :branches => begin
                rule = HerbCore.get_rule(node)
                @match solver.grammar.rules[rule] begin
                    :(branch * branches) => begin
                        # Get the relevant paths for the branch and branches
                        get_relevant_counts(
                            solver, c, child(path, 1), node.children[1], counts
                        )
                        get_relevant_counts(
                            solver, c, push!(path, 2), node.children[2], counts
                        )

                        return counts
                    end
                    :("") => return counts
                    x => throw("Unknown branches rule: $x")
                end
            end
            :branch => begin
                get_relevant_counts(solver, c, push!(path, 2), node.children[2], counts)
                return counts
            end
            x => begin
                atom_counts = c.rule_to_dict[get_rules(node)[1]]
                for (elem, count) in atom_counts
                    counts.atoms_fixed[elem] = get(counts.atoms_fixed, elem, 0) + count
                end
                return counts
            end
        end
    else
        @match node_type begin
            :ringbonds => begin
                return counts
            end
            :branches ||
            :chain ||
            :molecule_list ||
            :molecule => begin
                counts.holes += 1
                push!(counts.hole_paths, path)
                return counts
            end
            x => begin
                # return RelevantCounts(
                #     0, Dict{String, Int}(), 1, [path]
                # )
                atom_counts = c.rule_to_dict[get_rules(node)[1]]
                for (elem, count) in atom_counts
                    counts.atoms_fixed[elem] = get(counts.atoms_fixed, elem, 0) + count
                end
                return counts
            end
        end
    end
end

function HerbConstraints.propagate!(
        solver::Solver, constraint::LocalGenericBalancedReaction
)
    node = get_node_at_location(solver, constraint.path)
    if !isuniform(node)
        return nothing
    end

    left_path = child(constraint.path, 1)
    right_path = child(constraint.path, 2)

    left_node = get_node_at_location(solver, left_path)
    right_node = get_node_at_location(solver, right_path)
    
    # Get the relevant counts for the left side
    left_counts = get_relevant_counts(
        solver, constraint, left_path, left_node
    )

    # Only continue if the are no holes anymore on the left side (leftmost heuristic)
    if left_counts.holes > 0
        return
    end

    # Get the relevant counts for the right side
    right_counts = get_relevant_counts(
        solver, constraint, right_path, right_node
    )

    # Get the relevant atom counts
    left_atoms = left_counts.atoms_fixed
    right_atoms = right_counts.atoms_fixed
    all_keys = union(collect(keys(left_atoms)), collect(keys(right_atoms)))
    difference = Dict{String, Int}()
    for k in all_keys
        left_count = get(left_atoms, k, 0)
        right_count = get(right_atoms, k, 0)
        if left_count != right_count
            difference[k] = left_count - right_count
        end
    end

    # If there are any negative values in the difference (the left side needs additional atoms), the solver is infeasible
    if any(v -> v < 0, values(difference))
        HerbConstraints.set_infeasible!(solver)
        return
    end

    # If there are no holes on the right side and the right side needs additional atoms, the solver is infeasible
    if right_counts.holes == 0 && any(v -> v > 0, values(difference))
        HerbConstraints.set_infeasible!(solver)
        return
    end



    if length(keys(difference)) == 0
        for hole in right_counts.hole_paths
            node = get_node_at_location(solver, hole)
            c_remove_all_but!(solver, hole, constraint.hole_terminals, false)
        end
        return
    end

    if right_counts.holes == 1
        c_remove!(solver, right_counts.hole_paths[1], constraint.hole_terminals, false)
        return
    end


    # Gather the rules that make the solver infeasible
    # (Rules that make the atom counts on the right side more than the left side)
    impossible_rules = Vector{Int}()
    for (rule, count) in constraint.rule_to_dict
        for (k, rule_count) in count
            diff = get(difference, k, 0)
            if !(diff + rule_count > 0)
                push!(impossible_rules, rule)
                continue
            end
        end
    end
    
    for hole in right_counts.hole_paths
        c_remove!(solver, hole, impossible_rules, false)
    end
end
