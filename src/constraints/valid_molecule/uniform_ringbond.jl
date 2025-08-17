mutable struct localUniformRingbonds <: AbstractLocalConstraint
    path::Vector{Int}
    ringbond_paths::Vector{Tuple{Vector{Int}, Vector{Int}}}
    connected_groups::Vector{Vector{Vector{Int}}}
    grammar_data::GrammarData
end

import Base.==
function ==(a::localUniformRingbonds, b::localUniformRingbonds)
    return a.path == b.path && a.ringbonds == b.ringbonds
end

import Base.hash
function hash(a::localUniformRingbonds, h::UInt)
    return hash(a.path, h) + hash(a.ringbond_paths, h)
end

function HerbConstraints.shouldschedule(
        solver::Solver, constraint::localUniformRingbonds, path::Vector{Int}
)::Bool
    # Check if the update was in a ringbond
    for ringbond_path in [x[1] for x in constraint.ringbond_paths]
        if path[1:(end - 1)] == ringbond_path
            return true
        end
    end

    return false
end

function HerbConstraints.propagate!(solver::Solver, constraint::localUniformRingbonds)

    # TODO: make dynamic; TEMP: restrict max ringbonds digit to 9
    max_ringbond_grammar = 2
    # Constraint the max digit of the ringbonds
    max_ringbonds = length(constraint.ringbond_paths)÷2
    for (i, (ringbond_path, atom_path)) in enumerate(constraint.ringbond_paths)
        rule = digit_to_grammar(
            constraint.grammar_data, min(max_ringbonds, max_ringbond_grammar)
        )
        grammar_i = digit_to_grammar(constraint.grammar_data, min(i, max_ringbond_grammar))
        node = get_node_at_location(solver, push!(copy(ringbond_path), 2))
        if node isa StateHole
            c_remove_above!(
                solver, push!(copy(ringbond_path), 2), min(rule, grammar_i); fix_point = false
            )
        end
    end

    # Get the filled ringbond digits
    filled_ringbonds = Dict{Int, Tuple}()
    for (ringbond_path, atom_path) in constraint.ringbond_paths
        digit = get_node_at_location(solver, push!(copy(ringbond_path), 2))
        if isfilled(digit)
            rule = HerbCore.get_rule(digit)
            if haskey(filled_ringbonds, rule)
                if !isnothing(filled_ringbonds[rule][2])
                    HerbConstraints.set_infeasible!(solver)
                end

                (first_rb, _) = filled_ringbonds[rule]
                filled_ringbonds[rule] = (first_rb, ringbond_path)
            else
                filled_ringbonds[rule] = (ringbond_path, nothing)
            end
        end
    end

    # If the ringbond digit is the same, they both should have the same bond
    for (rule, pair) in collect(filled_ringbonds)
        if !isnothing(pair[2])
            # Get the bond of the first ringbond
            left_path = push!(copy(pair[1]), 1)
            left_node = get_node_at_location(solver, left_path)
            left_rules = get_rules(left_node)

            # Get the bond of the second ringbond
            right_path = push!(copy(pair[2]), 1)
            right_node = get_node_at_location(solver, right_path)
            right_rules = get_rules(right_node)

            # Make their domains equal
            c_remove_all_but!(solver, left_path, right_rules, false)
            c_remove_all_but!(solver, right_path, left_rules, false)
        end
    end

    # Check that a group should not have the same digits
    for group in constraint.connected_groups
        seen = Set{Int}()
        # println("Group: ", group)
        # Filter empty paths from the group
        group = filter(path -> !isempty(path), group)
        for ringbond_path in group
            # println("Ringbond path: ", ringbond_path)
            digit = get_node_at_location(solver, push!(copy(ringbond_path), 2))
            # println("Digit: ", digit)
            if isfilled(digit)
                rule = HerbCore.get_rule(digit)
                if rule in seen
                    HerbConstraints.set_infeasible!(solver)
                    # println("Infeasible: ", rule)
                end
                push!(seen, rule)
            end
        end
        # println("Seen: ", seen)
    end

    # Check that there are no two ringbonds between the same two atoms
    filled_atoms = Dict{Int, Tuple}()
    for (ringbond_path, atom_path) in constraint.ringbond_paths
        digit = get_node_at_location(solver, push!(copy(ringbond_path), 2))
        if isfilled(digit)
            rule = HerbCore.get_rule(digit)
            if haskey(filled_atoms, rule)
                if !isnothing(filled_atoms[rule][2])
                    HerbConstraints.set_infeasible!(solver)
                end

                (first_rb, _) = filled_atoms[rule]
                filled_atoms[rule] = (first_rb, atom_path)
            else
                filled_atoms[rule] = (atom_path, nothing)
            end
        end
    end
    # println("Filled atoms: ", filled_atoms)
    filled_atom_pairs = Dict{Tuple, Int}()
    for (rule, (atom1, atom2)) in collect(filled_atoms)
        if !isnothing(atom2)
            highest_atom = max(atom1, atom2)
            lowest_atom = min(atom1, atom2)
            if haskey(filled_atom_pairs, (highest_atom, lowest_atom))
                HerbConstraints.set_infeasible!(solver)
                # println("Infeasible: ", rule, " ", (highest_atom, lowest_atom), " ", filled_atom_pairs[(highest_atom, lowest_atom)])
            else
                filled_atom_pairs[(highest_atom, lowest_atom)] = rule
            end
        end
    end
    # println("Filled pairs: ", filled_atom_pairs)
end

function get_ringbond_paths(
        solver::UniformSolver,
        path::Vector{Int};
        forbidden_group::Vector{Vector{Int}} = Vector{Vector{Int}}(),
        atom_path::Vector{Int} = Vector{Int}()
)::Tuple{
        Vector{Tuple{Vector{Int}, Vector{Int}}}, Vector{Vector{Vector{Int}}}, Vector{Vector{Int}}
}

    # Get the node at the specified path
    node = get_node_at_location(solver, path)
    type = get_node_type(solver.grammar, node)

    @match type begin
        :molecule => return get_ringbond_paths(solver, push!(copy(path), 1))
        :chain => begin
            rule = solver.grammar.rules[HerbCore.get_rule(node)]
            @match rule begin
                :(SMILES_combine_chain(bond,
                    structure,
                    chain)) => begin
                    structure_ringbonds, structure_forbidden,
                    ringbond_group = get_ringbond_paths(
                        solver,
                        push!(copy(path), 2),
                        forbidden_group = forbidden_group
                    )
                    chain_ringbonds, chain_forbidden,
                    _ = get_ringbond_paths(
                        solver, push!(copy(path), 3), forbidden_group = ringbond_group
                    )
                    return vcat(chain_ringbonds, structure_ringbonds),
                    vcat(chain_forbidden, structure_forbidden),
                    []
                end
                :(structure * bond *
                  chain) => begin
                    structure_ringbonds, structure_forbidden,
                    ringbond_group = get_ringbond_paths(
                        solver,
                        push!(copy(path), 1),
                        forbidden_group = forbidden_group
                    )
                    chain_ringbonds, chain_forbidden,
                    _ = get_ringbond_paths(
                        solver, push!(copy(path), 3), forbidden_group = ringbond_group
                    )
                    return vcat(chain_ringbonds, structure_ringbonds),
                    vcat(chain_forbidden, structure_forbidden),
                    []
                end
                :(atom *
                  ringbonds) => begin
                    ringbonds_paths, _,
                    _ = get_ringbond_paths(
                        solver, push!(copy(path), 2), atom_path = push!(copy(path), 1)
                    )
                    forbidden_group = vcat(forbidden_group, [x[1] for x in ringbonds_paths])
                    ringbonds_group = copy([x[1] for x in ringbonds_paths])
                    return ringbonds_paths, [forbidden_group], ringbonds_group
                end
                _ => throw("Unknown chain rule: $rule")
            end
        end
        :structure => begin
            ringbond_paths, _,
            _ = get_ringbond_paths(
                solver, push!(copy(path), 2), atom_path = push!(copy(path), 1)
            )
            forbidden_group = vcat(forbidden_group, [x[1] for x in ringbond_paths])
            branches_ringbonds, branches_forbidden,
            _ = get_ringbond_paths(
                solver,
                push!(copy(path), 3),
                forbidden_group = [x[1] for x in ringbond_paths]
            )
            return vcat(ringbond_paths, branches_ringbonds),
            vcat([forbidden_group], branches_forbidden),
            [x[1] for x in ringbond_paths]
        end
        :ringbonds => begin
            rule = solver.grammar.rules[HerbCore.get_rule(node)]
            @match rule begin
                :("") => return [], [], []
                :(ringbond *
                  ringbonds) => begin
                    ringbond_path = push!(copy(path), 1)
                    ringbond_paths, _,
                    _ = get_ringbond_paths(
                        solver, push!(copy(path), 2), atom_path = atom_path
                    )
                    return vcat([(ringbond_path, atom_path)], ringbond_paths), [], []
                end
                _ => throw("Unknown ringbonds rule: $rule")
            end
        end
        :branches => begin
            rule = solver.grammar.rules[HerbCore.get_rule(node)]
            @match rule begin
                :("") => return [], [], []
                :(branch *
                  branches) => begin
                    branch_ringbonds, branch_forbidden,
                    _ = get_ringbond_paths(
                        solver,
                        push!(copy(path), 1),
                        forbidden_group = forbidden_group
                    )
                    branches_ringbonds, branches_forbidden,
                    _ = get_ringbond_paths(
                        solver,
                        push!(copy(path), 2),
                        forbidden_group = forbidden_group
                    )
                    return vcat(branch_ringbonds, branches_ringbonds),
                    vcat(branch_forbidden, branches_forbidden),
                    []
                end
                _ => throw("Unknown branches rule: $rule")
            end
        end
        :branch => begin
            return get_ringbond_paths(
                solver, push!(copy(path), 2), forbidden_group = forbidden_group
            )
        end
        _ => throw("Unknown node type: $type")
    end
end

function post_ringbond_constraints!(
        solver::Solver, path::Vector{Int}, grammar_data::GrammarData
)
    # Get the ringbonds underneath the current path
    ringbond_paths, incompatible_groups, _ = get_ringbond_paths(solver, path)

    # Check if there are an even number of ringbonds
    if length(ringbond_paths) % 2 != 0
        HerbConstraints.set_infeasible!(solver)
    end

    # Post the local constraint with the incompatible pairs
    HerbConstraints.post!(
        solver,
        localUniformRingbonds(path, ringbond_paths, incompatible_groups, grammar_data)
    )
end
