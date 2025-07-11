# TODO: check the place of this definition
# Define chemical properties
atom_valences = Dict(
    "[O]" => 2,
    "[H]" => 1,
    "[C]" => 4,
    "[N]" => 3
)
bond_orders = Dict(
    "-" => 1,
    "=" => 2,
    "≡" => 3,
    "≣" => 4
)
digits = Dict(
    "1" => 1,
    "2" => 2,
    "3" => 3,
    "4" => 4,
    "5" => 5,
    "6" => 6,
    "7" => 7,
    "8" => 8,
    "9" => 9
)


struct ValidSMILES <: AbstractGrammarConstraint
    grammar_data::GrammarData
end

function HerbCore.is_domain_valid(constraint::ValidSMILES, grammar::ContextSensitiveGrammar)
    # TODO: Implement the logic to check if the domain of the ValidSMILES constraint is valid
    return true
end

function HerbCore.update_rule_indices!(
    constraint::ValidSMILES, 
    n_rules::Integer,
    mapping::AbstractDict{<:Integer,<:Integer},
    constraints::Vector{<:AbstractConstraint}
)
    # TODO: Implement the logic to update the rule indices of the ValidSMILES constraint
    update_grammar_data!(constraint.grammar_data, mapping)
end

function HerbCore.update_rule_indices!(
    constraint::ValidSMILES, 
    n_rules::Integer
)
    # TODO: Implement the logic to update the rule indices of the ValidSMILES constraint
end

function HerbConstraints.on_new_node(solver::Solver, constraint::ValidSMILES, path::Vector{Int})
    if solver isa GenericSolver
        node = get_node_at_location(solver, path)
        type = get_node_type(solver.grammar, node)

        if type == :molecule
            constraint = GenericAtomConnections(path, constraint.grammar_data)
            HerbConstraints.post!(solver, constraint)
        end
    else
    
        node = get_node_at_location(solver, path)
        grammar = solver.grammar

        # Only proceed if node is filled and is a molecule rule
        if !isfilled(node) || grammar.types[HerbCore.get_rule(node)] != :molecule
            return
        end

        # Post the atom constraints
        post_atom_constraints!(solver, path, constraint.grammar_data)

        # Post the ringbond constraints
        post_ringbond_constraints!(solver, path, constraint.grammar_data)
    end
end





