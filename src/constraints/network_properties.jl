struct NetworkProperties <: AbstractGrammarConstraint 
    problem::ProblemDefinition
end


import Base.==
function ==(a::NetworkProperties, b::NetworkProperties)
    return a.problem == b.problem
end

import Base.hash
function hash(a::NetworkProperties, h::UInt)
    return hash(a.problem, h)
end



struct LocalNetworkProperties <: AbstractLocalConstraint
    path::Vector{Int}
    problem::ProblemDefinition
    contains_molecules::Dict{Molecule, StateInt}
    reaction_possible::Dict{AbstractLocalConstraint, Vector{Molecule}}
end

function LocalNetworkProperties(path, problem::ProblemDefinition, sm)
    contains_molecules = Dict{Molecule, StateInt}()
    for molecule in problem.known_molecules
        contains_molecules[molecule] = StateInt(sm, 0)
    end
    return LocalNetworkProperties(path, problem, contains_molecules, Dict{AbstractLocalConstraint, Vector{Molecule}}())
end


function get_reaction_paths(solver::Solver, path::Vector{Int})
    node = get_node_at_location(solver, path)
    type = get_node_type(solver.grammar, node)
    if type == :reaction
        return [path]
    end
    result = Vector{Vector{Int}}()
    for (i, child) in enumerate(node.children)
        next_path = push!(copy(path), i)
        result = vcat(result, get_reaction_paths(solver, next_path))
    end
    return result
end


function HerbConstraints.on_new_node(solver::Solver, constraint::NetworkProperties, path::Vector{Int})
    if solver isa GenericSolver
        return
    end
    
    node = get_node_at_location(solver, path)
    type = get_node_type(solver.grammar, node)

    if type == :network
        local_network_properties = LocalNetworkProperties(path, constraint.problem, solver.sm)

        reaction_paths = get_reaction_paths(solver, path)
        reaction_constraints = post_reaction_constraints!(solver, reaction_paths, local_network_properties)

        # HerbConstraints.post!(solver, local_network_properties)
    end
end



function HerbConstraints.propagate!(solver::Solver, constraint::LocalNetworkProperties)
    if solver isa GenericSolver
        return
    end

    # node = get_node_at_location(solver, constraint.path)
    # for molecule in keys(constraint.contains_molecules)
    #     if get_value(constraint.contains_molecules[molecule]) == 0
    #         # check if the molecule is in the reaction
    #         rule = findfirst(==(:($molecule)), solver.grammar.rules)
    #         @match _contains(node, rule) begin
    #             true => begin 
    #                 set_value!(constraint.contains_molecules[molecule], 1)
    #             end
    #             false => begin 
    #                 HerbConstraints.set_infeasible!(solver)
    #             end
    #             holes::Vector{AbstractHole} => begin
    #                 # @assert length(holes) > 0
    #                 # if length(holes) == 1
    #                 #     if isuniform(holes[1])
    #                 #         track!(solver, "LocalContains deduction")
    #                 #         path = vcat(c.path, get_path(node, holes[1]))
    #                 #         deactivate!(solver, c)
    #                 #         remove_all_but!(solver, path, c.rule)
    #                 #     else
    #                 #         # we cannot deduce anything yet, new holes can appear underneath this hole
    #                 #         # optimize this by checking if the target rule can appear as a child of the hole
    #                 #         track!(solver, "LocalContains softfail (non-uniform hole)")
    #                 #     end
    #                 # else
    #                 #     # multiple holes can be set to the target value, no deduction can be made as this point
    #                 #     # optimize by only repropagating if the number of holes involved is <= 2
    #                 #     track!(solver, "LocalContains softfail (>= 2 holes)")
    #                 # end
    #             end
    #         end
    #     end
    # end
end