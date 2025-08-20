struct BalancedReaction <: AbstractGrammarConstraint
    rule_to_dict::Dict{Int, Dict{String, Int}}
    hole_terminals::Vector{Int}  # Indices of terminal holes in the grammar
end

function BalancedReaction()
    return BalancedReaction(Dict{Int, Dict{String, Int}}(), Int[])
end

function HerbCore.is_domain_valid(
        constraint::BalancedReaction, grammar::ContextSensitiveGrammar
)
    for (i, rule) in enumerate(grammar.rules)
        if rule isa Molecule
            # Create a mapping from rule index to atom counts
            atom_counts = count_atoms(rule)
            constraint.rule_to_dict[i] = atom_counts
        end
    end

    hole_terminals = findall(grammar.isterminal)
    for terminal in hole_terminals
        if terminal ∉ constraint.hole_terminals
            push!(constraint.hole_terminals, terminal)
        end
    end

    # TODO: Check if the grammar is valid for balanced reactions
    return true
end

function HerbCore.update_rule_indices!(constraint::BalancedReaction, n_rules::Integer)
    return nothing
end

function HerbCore.update_rule_indices!(constraint::BalancedReaction, grammar::HerbGrammar.ContextSensitiveGrammar)
    for (i, rule) in enumerate(grammar.rules)
        if rule isa Molecule
            # Create a mapping from rule index to atom counts
            atom_counts = count_atoms(rule)
            constraint.rule_to_dict[i] = atom_counts
        end
    end

    hole_terminals = findall(grammar.isterminal)
    for terminal in hole_terminals
        if terminal ∉ constraint.hole_terminals
            push!(constraint.hole_terminals, terminal)
        end
    end
end

function HerbCore.update_rule_indices!(
        constraint::BalancedReaction,
        n_rules::Integer,
        mapping::AbstractDict{<:Integer, <:Integer},
        constraints::Vector{<:AbstractConstraint}
)
    return nothing
end

function HerbConstraints.on_new_node(
        solver::Solver, constraint::BalancedReaction, path::Vector{Int}
)
    node = get_node_at_location(solver, path)
    type = get_node_type(solver.grammar, node)

    if type != :reaction
        return
    end

    if solver isa UniformSolver
        HerbConstraints.post!(
            solver, LocalUniformBalancedReaction(solver, path)
        )
    else
        HerbConstraints.post!(
            solver, LocalGenericBalancedReaction(path, constraint.rule_to_dict, constraint.hole_terminals)
        )
    end
end

function is_valid(candidate::Reaction, constraint::BalancedReaction)
    input_counts = Dict{String, Int}()
    output_counts = Dict{String, Int}()

    for (num, molecule) in candidate.inputs
        for (atom, count) in count_atoms(molecule)
            input_counts[atom] = get(input_counts, atom, 0) + count * num
        end
    end

    for (num, molecule) in candidate.outputs
        for (atom, count) in count_atoms(molecule)
            output_counts[atom] = get(output_counts, atom, 0) + count * num
        end
    end

    if input_counts != output_counts
        return false
    end

    return true
end
