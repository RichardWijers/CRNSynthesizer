struct ContainsMolecules <: AbstractGrammarConstraint
    rules::Vector{Tuple{Int, ReactionPosition}}
end

struct LocalGenericContainsMolecules <: AbstractLocalConstraint
	path::Vector{Int}
    rules::Vector{Tuple{Int, ReactionPosition}}
end

struct LocalUniformContainsMolecules <: AbstractLocalConstraint
    path::Vector{Int}
    rules::Vector{Tuple{Int, ReactionPosition}}
    required_paths::Vector{Vector{Int}}
    required_path_groups::Vector{Vector{Vector{Int}}}
    possibilities::Vector{Vector{Int}}
end

function get_required_molecules_paths(solver::Solver, path::Vector{Int}, result::Vector{Vector{Vector{Int}}} = Vector{Vector{Vector{Int}}}())
    function get_paths(solver, path::Vector{Int}, result::Vector{Vector{Int}} = Vector{Vector{Int}}())
        node = get_node_at_location(solver, path)
        type = get_node_type(solver.grammar, node)
        if type == :required_molecule
            return [path]
        end
        for (i, child) in enumerate(node.children)
            next_path = push!(copy(path), i)
            result = vcat(result, get_paths(solver, next_path))
        end
        return result
    end

    node = get_node_at_location(solver, path)
    type = get_node_type(solver.grammar, node)
    if type == :reaction
        left_paths = get_paths(solver, push!(copy(path), 1))
        right_paths = get_paths(solver, push!(copy(path), 2))
        push!(result, left_paths)
        push!(result, right_paths)
    end
    for (i, child) in enumerate(node.children)
        next_path = push!(copy(path), i)
        get_required_molecules_paths(solver, next_path, result)
    end
    return result
end


# function search_required_molecules(solver::Solver, partial::Vector{Int}, index::Int, left_positions::Vector{Int}, right_positions::Vector{Int}, ordered::Vector{Tuple{Int,Int}}, required_left::Set{Int}, required_right::Set{Int}, required::Set{Int}, flattened::Vector{Vector{Int}}, valid_combinations::Vector{Vector{Int}})
#     if index > length(flattened)
#         # Check if this combination satisfies all constraints
#         left_values = partial[left_positions]
#         right_values = partial[right_positions]
        
#         if issubset(required, partial) && 
#             issubset(required_left, left_values) &&
#             issubset(required_right, right_values)
#             push!(valid_combinations, copy(partial))
#         end
#         return
#     end
    
#     for value in get_rules(get_node_at_location(solver, flattened[index]))
#         # Check if adding this value violates any ordering constraint
#         valid = true
#         for pair in ordered
#             x, y = pair[1], pair[2]
#             if y == index && x < index && partial[x] > value
#                 valid = false
#                 break
#             elseif x == index && y < index && value > partial[y]
#                 valid = false
#                 break
#             end
#         end
        
#         if valid
#             partial[index] = value
#             search_required_molecules(solver, partial, index + 1, left_positions, right_positions, ordered, required_left, required_right, required, flattened, valid_combinations)
#         end
#     end
# end

# function get_required_molecules_possibilities2(solver::Solver, paths::Vector{Vector{Vector{Int}}}, flattened::Vector{Vector{Int}}, rules::Vector{Tuple{Int, ReactionPosition}})::Vector{Vector{Int}}
#     required = Set{Int}()
#     required_left = Set{Int}()
#     required_right = Set{Int}()
#     for (rule, position) in rules
#         if position == INPUT
#             push!(required_left, rule)
#         elseif position == OUTPUT
#             push!(required_right, rule)
#         end
#         push!(required, rule)
#     end

#     left_positions = Vector{Int}()
#     right_positions = Vector{Int}()
#     ordered = Vector{Tuple{Int,Int}}()
#     position = 1
#     for (i, path_group) in enumerate(paths)
#         if isodd(i)
#             append!(left_positions, position:(position+length(path_group)-1))
#         else
#             append!(right_positions, position:(position+length(path_group)-1))
#         end
#         for j in 0:(length(path_group)-2)
#             push!(ordered, (position+j, position+j+1))
#         end
#         position += length(path_group)
#     end
    
#     valid_combinations = Vector{Vector{Int}}()
    
    

#     initial_partial::Vector{Int} = zeros(Int, length(flattened))
#     search_required_molecules(solver, initial_partial, 1, left_positions, right_positions, ordered, required_left, required_right, required, flattened, valid_combinations)

#     return valid_combinations
# end


function get_required_molecules_possibilities(solver::Solver, paths::Vector{Vector{Vector{Int}}}, flattened::Vector{Vector{Int}}, rules::Vector{Tuple{Int, ReactionPosition}})::Vector{Vector{Int}}
    # Initialize with empty sequence
    current_options = [Int[]]
    
    # Calculate reaction boundaries for quick lookup
    reaction_boundaries = [0]
    for path_group in paths
        push!(reaction_boundaries, reaction_boundaries[end] + length(path_group))
    end
    
    # Process each molecule position
    for (reaction_idx, path_group) in enumerate(paths)
        for path in path_group
            node = get_node_at_location(solver, path)
            rules::Vector{Int} = get_rules(node)
            
            new_options = Vector{Vector{Int}}()
            
            for option in current_options
                # Check if this is the start of a new reaction
                is_new_reaction = length(option) == reaction_boundaries[reaction_idx]
                # is_new_reaction == false
                
                for rule in rules
                    # Enforce non-decreasing order within a reaction
                    if is_new_reaction || rule >= option[end]
                        push!(new_options, vcat(option, rule))
                    end
                end
            end
            
            current_options = new_options
        end
    end
    
    return current_options
end

function HerbConstraints.on_new_node(solver::Solver, c::ContainsMolecules, path::Vector{Int})
    if length(path) == 0

        if solver isa GenericSolver
            HerbConstraints.post!(solver, LocalGenericContainsMolecules(path, c.rules))
        else
            required_paths = get_required_molecules_paths(solver, path)
            flattened_required_paths = Vector{Vector{Int}}()
            for path_group in required_paths
                for path in path_group
                    push!(flattened_required_paths, path)
                end
            end

            # First pass for caching the required paths
            rules = [x[1] for x in c.rules]
            possibilities::Vector{Vector{Int}} = []
            # possibilities::Vector{Vector{Int}} = get_required_molecules_possibilities(solver, required_paths, flattened_required_paths, c.rules)
            # filter!(x -> issubset(rules, x), possibilities)

            # left_positions::Vector{Int} = Int[]
            # right_positions::Vector{Int} = Int[]
            # left = true
            # position = 1
            # for path_group in required_paths
            #     if left
            #         left_positions = vcat(left_positions, position:(position+length(path_group)-1))
            #     else 
            #         right_positions = vcat(right_positions, position:(position+length(path_group)-1))
            #     end
            #     position += length(path_group)
            #     left = !left
            # end

            # for (rule, position) in c.rules
            #     if position == INPUT
            #         filter!(x -> rule in x[left_positions], possibilities)
            #     elseif position == OUTPUT
            #         filter!(x -> rule in x[right_positions], possibilities)
            #     end
            # end

            # if length(possibilities) == 0
            #     HerbConstraints.set_infeasible!(solver)
            # end

            HerbConstraints.post!(solver, LocalUniformContainsMolecules(path, c.rules, flattened_required_paths, required_paths, possibilities))
        end
    end
end


mutable struct RequiredPlaces
    input_holes::Vector{AbstractHole}
    output_holes::Vector{AbstractHole}
    input_count::Int
    output_count::Int
    non_uniform::Bool
    satisfied::Bool

    function RequiredPlaces(;
        input_holes::Vector{AbstractHole}=Vector{AbstractHole}(), 
        output_holes::Vector{AbstractHole}=Vector{AbstractHole}(), 
        input_count::Int=0, 
        output_count::Int=0, 
        non_uniform::Bool=false,
        satisfied::Bool=false
    )
        new(input_holes, output_holes, input_count, output_count, non_uniform, satisfied)
    end
end

function HerbConstraints.propagate!(solver::Solver, c::LocalGenericContainsMolecules)
    node = get_node_at_location(solver, c.path)
    input_rules = [rule for (rule, position) in c.rules if position == INPUT]
    output_rules = [rule for (rule, position) in c.rules if position == OUTPUT]
    unknown_rules = [rule for (rule, position) in c.rules if position == UNKNOWN]


    found_places = _find_required_holes(solver, node, input_rules, output_rules, unknown_rules)


    terminal_node = findfirst(==(:(Vector{Molecule}())), solver.grammar.rules)



    if found_places.satisfied
        HerbConstraints.deactivate!(solver, c)
        return
    end

    if found_places.non_uniform
        return
    end

    # println(node)
    # println("Found places: ", found_places)
    # println("Terminal node: ", terminal_node)
    # println("Input rules: ", input_rules)
    # println("Output rules: ", output_rules)
    # println("Unknown rules: ", unknown_rules)

    if length(input_rules) > found_places.input_count 
        if length(found_places.input_holes) == 0
            HerbConstraints.set_infeasible!(solver)
        elseif length(found_places.input_holes) == 1
            path = get_path(node, found_places.input_holes[1])
            remove!(solver, path, terminal_node)
        end
    end

    if length(output_rules) > found_places.output_count 
        if length(found_places.output_holes) == 0
            HerbConstraints.set_infeasible!(solver)
        elseif length(found_places.output_holes) == 1
            path = get_path(node, found_places.output_holes[1])
            remove!(solver, path, terminal_node)
        end
    end

    if length(found_places.input_holes) == 0 && length(found_places.output_holes) == 0
        HerbConstraints.set_infeasible!(solver)
    end

    if length(found_places.input_holes) == 1 && length(found_places.output_holes) == 0
        path = get_path(node, found_places.input_holes[1])
        remove!(solver, path, terminal_node)
    end
    if length(found_places.input_holes) == 0 && length(found_places.output_holes) == 1
        path = get_path(node, found_places.output_holes[1])
        remove!(solver, path, terminal_node)
    end

end


function HerbConstraints.shouldschedule(solver::Solver, c::LocalUniformContainsMolecules, path::Vector{Int})
    if path in c.required_paths
        return true
    end
    return false
end


function HerbConstraints.propagate!(solver::Solver, c::LocalUniformContainsMolecules)
    # Propagation when the list of possibilities is filled
    # possibilities = copy(c.possibilities)
    # for (i, path) in enumerate(c.required_paths)
    #     node = get_node_at_location(solver, path)
    #     rules = get_rules(node)
    #     filter!(x -> x[i] in rules, possibilities)
    # end

    # for (i, path) in enumerate(c.required_paths)
    #     possible = [x[i] for x in possibilities]
    #     c_remove_all_but!(solver, path, possible, false)
    #     if !isfeasible(solver)
    #         return
    #     end
    # end



    required_input_rules = [rule for (rule, position) in c.rules if position == INPUT]
    required_output_rules = [rule for (rule, position) in c.rules if position == OUTPUT]
    required_unknown_rules = [rule for (rule, position) in c.rules if position == UNKNOWN]

    input_rules = Set()
    input_holes = 0
    output_rules = Set()
    output_holes = 0
    all_rules = Set()
    for (group, paths) in enumerate(c.required_path_groups)
        if isodd(group)
            for path in paths
                node = get_node_at_location(solver, path)
                if isfilled(node)
                    rule = get_rule(node)
                    if rule in required_input_rules
                        push!(input_rules, get_rule(node))
                    end
                    push!(all_rules, get_rule(node))
                else
                    input_holes += 1
                end
            end
        else
            for path in paths
                node = get_node_at_location(solver, path)
                if isfilled(node)
                    rule = get_rule(node)
                    if rule in required_output_rules
                        push!(output_rules, get_rule(node))
                    end
                    push!(all_rules, get_rule(node))
                else
                    output_holes += 1
                end
            end
        end
    end

    # println("Input rules: ", input_rules, " (", length(input_rules), ")")
    # println("Output rules: ", output_rules, " (", length(output_rules), ")")
    # println("Input holes: ", input_holes)
    # println("Output holes: ", output_holes)
    # println("All rules: ", all_rules, " (", length(all_rules), ")")
    # println("Required input rules: ", required_input_rules, " (", length(required_input_rules), ")")
    # println("Required output rules: ", required_output_rules, " (", length(required_output_rules), ")")
    # println("Required unknown rules: ", required_unknown_rules, " (", length(required_unknown_rules), ")")
    # println()

    unknown_rules = filter(x -> !(x in required_input_rules || x in required_output_rules), required_unknown_rules)
    if length(required_input_rules) > length(input_rules) + input_holes || 
        length(required_output_rules) > length(output_rules) + output_holes
        length(required_input_rules) + length(required_output_rules) + length(unknown_rules) > length(input_rules) + length(output_rules) + input_holes + output_holes
        HerbConstraints.set_infeasible!(solver)
        return
    end
end

function _find_required_holes(
        solver::Solver,
        node::Union{AbstractHole, AbstractUniformHole, AbstractRuleNode},
        input_rules::Vector{Int},
        output_rules::Vector{Int},
        unknown_rules::Vector{Int},
        found_places::RequiredPlaces=RequiredPlaces(),
        current_position::ReactionPosition = UNKNOWN
    )
    node_type = get_node_type(solver.grammar, node)
    required_rule = findfirst(==(:(vcat([required_molecule], molecule_list))), solver.grammar.rules)

    @match node_type begin
        :molecule_list => begin
            r = get_rules(node)
            if required_rule in r
                if !isuniform(node)
                    if current_position == INPUT
                        push!(found_places.input_holes, node)
                    elseif current_position == OUTPUT
                        push!(found_places.output_holes, node)
                    else
                        throw(ArgumentError("Position $current_position is not valid"))
                    end
                end
            end
            for child::Union{AbstractHole, AbstractUniformHole, AbstractRuleNode} in get_children(node)
                _find_required_holes(solver, child, input_rules, output_rules, unknown_rules, found_places, current_position)
            end
        end

        :network || :reaction_list  => begin 

            if !isuniform(node)
                found_places.non_uniform = true
                return found_places
            end

            for child::Union{AbstractHole, AbstractUniformHole, AbstractRuleNode} in get_children(node)
                _find_required_holes(solver, child, input_rules, output_rules, unknown_rules, found_places, current_position)
            end
        end

        :reaction => begin
            if !isuniform(node)
                found_places.non_uniform = true
                return found_places
            end

            children::Vector{Union{AbstractHole, AbstractUniformHole, AbstractRuleNode}} = get_children(node)

            _find_required_holes(solver, children[1], input_rules, output_rules, unknown_rules, found_places, INPUT)
            _find_required_holes(solver, children[2], input_rules, output_rules, unknown_rules, found_places, OUTPUT)
        end

        :molecule => begin end

        :required_molecule => begin
            if current_position == INPUT
                found_places.input_count += 1
            elseif current_position == OUTPUT 
                found_places.output_count += 1
            else
                throw(ArgumentError("Position $current_position is not valid"))
            end

            if found_places.input_count >= length(input_rules) && 
                found_places.output_count >= length(output_rules) && 
                found_places.input_count + found_places.output_count >= length(input_rules) + length(output_rules) + length(unknown_rules)
                found_places.satisfied = true
            end
        end

        _ => throw(ArgumentError("Unknown rule type: $node_type"))
    end

    return found_places
end