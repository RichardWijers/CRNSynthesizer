struct BalancedReaction <: AbstractGrammarConstraint 
    complete_grammar::Bool
    problem::ProblemDefinition

    function BalancedReaction(; complete_grammar::Bool=true, problem::ProblemDefinition=ProblemDefinition())
        return new(complete_grammar, problem)
    end
end


struct AtomCounts
    atom_counts::Dict{String, Int}
    hash::UInt
end
function AtomCounts(atom_counts::Dict{String, Int})
    return AtomCounts(atom_counts, hash(atom_counts))
end

function Base.show(io::IO, a::AtomCounts)
    atom_counts = a.atom_counts
    atom_counts_str = join([string(k) * (v == 1 ? "" : string(v)) for (k, v) in atom_counts], " + ")
    print(io, atom_counts_str)
end

import Base.==
function ==(a::AtomCounts, b::AtomCounts)
    return a.hash == b.hash
end

import Base.hash
function hash(a::AtomCounts, h::UInt)
    return a.hash
end


struct ReactionPossibility
    rules::Vector{Int}
    atom_counts::AtomCounts
    hash::UInt
end
function ReactionPossibility(rules::Vector{Int}, atom_counts::Dict{String, Int})
    atom_counts = AtomCounts(atom_counts)
    return ReactionPossibility(rules, atom_counts, hash(rules, hash(atom_counts)))
end
function ReactionPossibility(rules::Vector{Int}, atom_counts::AtomCounts)
    return ReactionPossibility(rules, atom_counts, hash(rules, hash(atom_counts)))
end

function Base.show(io::IO, a::ReactionPossibility)
    rules_str = join([string(r) for r in a.rules], ", ")
    atom_counts_str = string(a.atom_counts)
    print(io, "Rules: [", rules_str, "], Atom Counts: ", atom_counts_str)
end

import Base.==
function ==(a::ReactionPossibility, b::ReactionPossibility)
    return a.hash == b.hash
end

import Base.hash
function hash(a::ReactionPossibility, h::UInt)
    return a.hash
end


struct LocalBalancedReaction <: AbstractLocalConstraint 
    path::Vector{Int}
    input_molecule_paths::Vector{Vector{Int}}
    output_molecule_paths::Vector{Vector{Int}}
    rule_to_atoms::Dict{Int, Dict{Any, Any}}
    left_possibilities::Vector{ReactionPossibility}
    right_possibilities::Vector{ReactionPossibility}
    NetworkProperties::LocalNetworkProperties
    last_reaction::Bool
end
function LocalBalancedReaction(path::Vector{Int}, input_molecule_paths::Vector{Vector{Int}}, output_molecule_paths::Vector{Vector{Int}}, rule_to_atoms::Dict{Int, Dict{Any, Any}}, left_possibilities::Vector{ReactionPossibility}, right_possibilities::Vector{ReactionPossibility}, last_reaction::Bool)
    local_network_properties = LocalNetworkProperties(Vector{Int}(), ProblemDefinition(), Dict{Molecule, StateInt}(), Dict{AbstractLocalConstraint, Vector{Molecule}}())
    return LocalBalancedReaction(path, input_molecule_paths, output_molecule_paths, rule_to_atoms, left_possibilities, right_possibilities, local_network_properties, last_reaction)
end

import Base.==
function ==(a::LocalBalancedReaction, b::LocalBalancedReaction)
    return a.path == b.path && a.input_molecule_paths == b.input_molecule_paths && a.output_molecule_paths == b.output_molecule_paths && a.last_reaction == b.last_reaction
end

import Base.hash
function hash(a::LocalBalancedReaction, h::UInt)
    return hash(a.path, h) + hash(a.input_molecule_paths, h) + hash(a.output_molecule_paths, h) + hash(a.last_reaction, h)
end





function get_atom_paths(solver::Solver, path::Vector{Int})
    node = get_node_at_location(solver, path)
    type = get_node_type(solver.grammar, node)
    if type == :atom
        return [path]
    end
    result = Vector{Vector{Int}}()
    for (i, child) in enumerate(node.children)
        next_path = push!(copy(path), i)
        result = vcat(result, get_atom_paths(solver, next_path))
    end
    return result
end

function get_molecule_paths(solver::Solver, path::Vector{Int})
    node = get_node_at_location(solver, path)
    type = get_node_type(solver.grammar, node)
    if type == :molecule
        return [path]
    end
    if type == :required_molecule
        return [path]
    end
    result = Vector{Vector{Int}}()
    for (i, child) in enumerate(node.children)
        next_path = push!(copy(path), i)
        result = vcat(result, get_molecule_paths(solver, next_path))
    end
    return result
end

function post_reaction_constraints!(solver::Solver, reaction_paths::Vector{Vector{Int}}, local_network_properties::LocalNetworkProperties)

    reaction_constraints = []
    for (i, path) in enumerate(reaction_paths)
        input_paths = get_molecule_paths(solver, push!(copy(path), 1))
        output_paths = get_molecule_paths(solver, push!(copy(path), 2))
        
        # updated_paths = Vector{Vector{Int}}()
        # HerbConstraints.post!(solver, LocalBalancedReaction(path, input_paths, output_paths, updated_paths, possible_left, possible_right))
        # HerbConstraints.post!(solver, LocalBalancedReaction(path, input_paths, output_paths))

        rule_to_atoms = Dict{Int, Dict{Any, Any}}()
        molecule_rules = findall(solver.grammar.domains[:molecule])

        for rule in molecule_rules
            rule_to_atoms[rule] = count_atoms(solver.grammar.rules[rule])
        end

        possible_left = get_possibilities(solver, input_paths, rule_counts=rule_to_atoms)
        possible_right = get_possibilities(solver, output_paths, rule_counts=rule_to_atoms)

        left_length = length(possible_left)
        right_length = length(possible_right)

        # Do a pre pass with the intersection
        l = [x.atom_counts for x in possible_left]
        r = [x.atom_counts for x in possible_right]
        intersection = intersect(l, r)
        filter!(x -> x.atom_counts in intersection, possible_left)
        filter!(x -> x.atom_counts in intersection, possible_right)

        if left_length == 0 || right_length == 0
            HerbConstraints.set_infeasible!(solver)
            return
        end


        balanced_reaction = LocalBalancedReaction(
            path, 
            input_paths, 
            output_paths, 
            rule_to_atoms, 
            possible_left, 
            possible_right, 
            local_network_properties,
            i == length(reaction_paths)
        )
        HerbConstraints.post!(solver, balanced_reaction)
        push!(reaction_constraints, balanced_reaction)
    end

    return reaction_constraints
end

function HerbConstraints.on_new_node(solver::Solver, constraint::BalancedReaction, path::Vector{Int})
    node = get_node_at_location(solver, path)
    type = get_node_type(solver.grammar, node)

    if constraint.complete_grammar
        if solver isa GenericSolver
            
            if type == :reaction
                # TODO: check the starting vectors for the input and output paths
                HerbConstraints.post!(solver, LocalGenericBalancedReaction(path))
            end
        else
            # Uniform solver
            node = get_node_at_location(solver, path)
            type = get_node_type(solver.grammar, node)
            if type == :reaction
                input_paths = get_atom_paths(solver, push!(copy(path), 1))
                output_paths = get_atom_paths(solver, push!(copy(path), 2))
                HerbConstraints.post!(
                    solver, LocalUniformBalancedReaction(path, input_paths, output_paths))
            end
        end
    else
        if solver isa UniformSolver && type == :reaction
            input_paths = get_molecule_paths(solver, push!(copy(path), 1))
            output_paths = get_molecule_paths(solver, push!(copy(path), 2))
            
            # updated_paths = Vector{Vector{Int}}()
            # HerbConstraints.post!(solver, LocalBalancedReaction(path, input_paths, output_paths, updated_paths, possible_left, possible_right))
            # HerbConstraints.post!(solver, LocalBalancedReaction(path, input_paths, output_paths))

            rule_to_atoms = Dict{Int, Dict{Any, Any}}()
            molecule_rules = findall(solver.grammar.domains[:molecule])

            for rule in molecule_rules
                rule_to_atoms[rule] = count_atoms(solver.grammar.rules[rule])
            end

            possible_left = get_possibilities(solver, input_paths, rule_counts=rule_to_atoms)
            possible_right = get_possibilities(solver, output_paths, rule_counts=rule_to_atoms)

            left_length = length(possible_left)
            right_length = length(possible_right)

            # Do a pre pass with the intersection
            l = [x.atom_counts for x in possible_left]
            r = [x.atom_counts for x in possible_right]
            intersection = intersect(l, r)
            filter!(x -> x.atom_counts in intersection, possible_left)
            filter!(x -> x.atom_counts in intersection, possible_right)

            if left_length == 0 || right_length == 0
                HerbConstraints.set_infeasible!(solver)
                return
            end

            HerbConstraints.post!(solver, LocalBalancedReaction(path, input_paths, output_paths, rule_to_atoms, possible_left, possible_right, false))
        end
    end
end


function HerbConstraints.shouldschedule(solver::Solver, constraint::LocalBalancedReaction, path::Vector{Int})
    if path in constraint.input_molecule_paths || path in constraint.output_molecule_paths
        # push!(constraint.updated_paths, path)
        return true
    end

    return false
end


function count_atoms(solver::Solver, paths::Vector{Vector{Int}})::Dict{String, Int}
    result = Dict{String, Int}()
    for path in paths
        node = get_node_at_location(solver, path)
        type = get_node_type(solver.grammar, node)

        if type == :molecule 
            if isfilled(node)
                rule = solver.grammar.rules[HerbCore.get_rule(node)]
                result = mergewith(+, result, count_atoms(rule))
            else 
                return nothing
            end
        end
    end

    return result
end

function get_possibilities(solver::Solver, paths::Vector{Vector{Int}}; rule_counts::Dict{Int, Dict{Any, Any}}=Dict{Int, Dict{Any, Any}}())

    # Test
    # indices: H => 1, O => 2, N => 3, C => 4
    atom_indices::Dict{String, Int} = Dict("H" => 1, "O" => 2, "N" => 3, "C" => 4)

    # Initialize with a single empty option
    current_options = [(Int[0,0,0,0], Int[])]
    
    # Process each path
    for path in paths
        node = get_node_at_location(solver, path)
        rules = get_rules(node)
        
        new_options = Vector{Tuple{Vector{Int}, Vector{Int}}}()
        
        for (atom_counts, rules_used) in current_options
            for rule_ind in rules
                # Enforce non-decreasing order
                if isempty(rules_used) || rule_ind >= rules_used[end]
                    # Get or calculate atom counts for this rule
                    if !haskey(rule_counts, rule_ind)
                        count = count_atoms(solver.grammar.rules[rule_ind])
                        rule_counts[rule_ind] = count
                    else
                        count = rule_counts[rule_ind]
                    end
                    
                    # Merge atom counts
                    merged_count::Vector{Int} = copy(atom_counts)
                    for (atom, count) in count
                        merged_count[atom_indices[atom]] += count
                    end

                    merged_rules::Vector{Int} = vcat(rules_used, rule_ind)
                    
                    push!(new_options, (merged_count, merged_rules))
                end
            end
        end
        
        current_options = new_options
    end
    
    # Convert to final result format
    results = ReactionPossibility[]
    
    for (counts, rules) in current_options
        atom_counts = Dict{String, Int}()
        for (i, count) in enumerate(counts)
            if count > 0
                if i == 1
                    atom_counts["H"] = count
                elseif i == 2
                    atom_counts["O"] = count
                elseif i == 3
                    atom_counts["N"] = count
                elseif i == 4
                    atom_counts["C"] = count
                else
                    error("Unexpected atom index: $i")
                end
            end
        end
        push!(results, ReactionPossibility(rules, atom_counts))
    end

    return results
end


function HerbConstraints.propagate!(solver::Solver, constraint::LocalBalancedReaction)
    if solver isa GenericSolver
        return
    end

    possible_left = copy(constraint.left_possibilities)
    for (index, path) in enumerate(constraint.input_molecule_paths)
        node = get_node_at_location(solver, path)
        rules = Set{Int}(get_rules(node))
        filter!(x -> x.rules[index] in rules, possible_left)
    end

    possible_right = copy(constraint.right_possibilities)
    for (index, path) in enumerate(constraint.output_molecule_paths)
        node = get_node_at_location(solver, path)
        rules = Set{Int}(get_rules(node))
        filter!(x -> x.rules[index] in rules, possible_right)
    end

    
    l = Set([x.atom_counts for x in possible_left])
    r = Set([x.atom_counts for x in possible_right])
    intersection = intersect(l, r)

    if isempty(intersection)
        HerbConstraints.set_infeasible!(solver)
        return
    end

    filter!(x -> x.atom_counts in intersection, possible_left)
    filter!(x -> x.atom_counts in intersection, possible_right)

    if isempty(possible_left) || isempty(possible_right)
        HerbConstraints.set_infeasible!(solver)
        return
    end
    
    fixed_left = Int[]
    for (i, path) in enumerate(constraint.input_molecule_paths)
        possibilities = [x.rules[i] for x in possible_left]
        # println("Possibilities for input ", i, ": ", possibilities)
        c_remove_all_but!(solver, path, possibilities, false)
        if !isfeasible(solver)
            return
        end

        node = get_node_at_location(solver, path)
        rules = get_rules(node)
        if length(possibilities) == 1 || length(rules) == 1
            push!(fixed_left, possibilities[1])
        end
    end

    fixed_right = Int[]
    for (i, path) in enumerate(constraint.output_molecule_paths)
        possibilities = [x.rules[i] for x in possible_right]
        # println("Possibilities for output ", i, ": ", possibilities)
        c_remove_all_but!(solver, path, possibilities, false)

        node = get_node_at_location(solver, path)
        rules = get_rules(node)
        if length(possibilities) == 1 || length(rules) == 1
            push!(fixed_right, possibilities[1])
        end
    end

    # any item that is fixed in the right must not be on the left
    for path in constraint.input_molecule_paths
        c_remove!(solver, path, fixed_right, false)
    end

    # any item that is fixed in the left must not be on the right
    for path in constraint.output_molecule_paths
        c_remove!(solver, path, fixed_left, false)
    end

    if (length(possible_left) == 1)
        for (i, path) in enumerate(constraint.input_molecule_paths)
            c_remove_all_but!(solver, path, [possible_left[1].rules[i]], false)
        end
    end

    if (length(possible_right) == 1)
        for (i, path) in enumerate(constraint.output_molecule_paths)
            c_remove_all_but!(solver, path, [possible_right[1].rules[i]], false)
        end
    end

    if !isfeasible(solver)
        return
    end



    required_molecules = constraint.NetworkProperties.problem.known_molecules
    fixed_left = map(x -> solver.grammar.rules[x], fixed_left)
    fixed_right = map(x -> solver.grammar.rules[x], fixed_right)

    # fixed = vcat(fixed_left, fixed_right)
    # for f in fixed
    #     if f in required_molecules && get_value(constraint.NetworkProperties.contains_molecules[f]) != 1
    #         set_value!(constraint.NetworkProperties.contains_molecules[f], 1)
    #     end
    # end
end