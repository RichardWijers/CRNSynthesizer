struct ContainsReactions <: AbstractGrammarConstraint
    required_molecules::Dict{RequiredMolecule, Vector{Int}}
end

function HerbCore.is_domain_valid(constraint::ContainsReactions, grammar::ContextSensitiveGrammar)
    # TODO: Implement the logic to check if the domain of the ContainsReactions constraint is valid
    return true
end

struct LocalContainsReactions <: AbstractLocalConstraint
    path::Vector{Int}
    required_molecules::Dict{RequiredMolecule, Vector{Int}}
end

function HerbConstraints.on_new_node(solver::Solver, constraint::ContainsReactions, path::Vector{Int})
    if length(path) == 0 && solver isa UniformSolver

        # for (required_molecule, indices) in constraint.required_molecules
        #     println(indices)
        # end

        HerbConstraints.post!(solver, LocalContainsReactions(path, constraint.required_molecules))
    end
end


function HerbConstraints.propagate!(solver::Solver, constraint::LocalContainsReactions)



    # Check if there are holes
    found = falses(length(constraint.required_molecules))
    found, holes = check_contains_reactions(solver, constraint, push!(copy(constraint.path), 1), found)

    # Check if the constraint is satisfied
    if all(found)
        HerbConstraints.deactivate!(solver, constraint)
    elseif holes == 0
        HerbConstraints.set_infeasible!(solver)
        return
    elseif holes == 1
        # Get all the missing required molecule indices
        required_indices = nothing
        for (found_ind, (required_molecule, indices)) in enumerate(constraint.required_molecules)
            if !found[found_ind]
                if isnothing(required_indices)
                    required_indices = Set(indices)
                else
                    required_indices = intersect(required_indices, Set(indices))
                end
            end
        end
        if length(required_indices) == 0
            # If there are no matching indices, the constraint is infeasible
            HerbConstraints.set_infeasible!(solver)
            return
        end
    end

    # println("Found: ", found)


end


function check_contains_reactions(solver::Solver, constraint::LocalContainsReactions, path::Vector{Int}, found, holes=0)

    node = get_node_at_location(solver, path)
    rule = solver.grammar.rules[get_rule(node)]

    @match rule begin

        :(vcat(reaction, reaction_list)) => begin
            node = get_node_at_location(solver, push!(copy(path), 1))
            if !isfilled(node)
                holes += 1
            else
                for (found_ind, (required_molecule, indices)) in enumerate(constraint.required_molecules)
                    if get_rule(node) in indices
                        found[found_ind] = true
                    end
                end
            end
            return check_contains_reactions(solver, constraint, push!(copy(path), 2), found, holes)
        end

        :([reaction]) => begin
            node = get_node_at_location(solver, push!(copy(path), 1))
            if !isfilled(node)
                holes += 1
            else
                for (found_ind, (required_molecule, indices)) in enumerate(constraint.required_molecules)
                    if get_rule(node) in indices
                        found[found_ind] = true
                    end
                end
            end
            return found, holes
        end

        :(Vector{Reaction}()) => begin
            return found, holes
        end

        _ => throw("Unexpected rule: $rule")
    end
end