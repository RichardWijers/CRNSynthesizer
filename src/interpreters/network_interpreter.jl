function interpret_network(
        program::AbstractRuleNode, grammar::AbstractGrammar
)::ReactionNetwork
    reactions::Vector{Reaction} = interpret_network_reactions(program.children[1], grammar)
    return ReactionNetwork(reactions)
end

function interpret_network_reactions(
        program::AbstractRuleNode, grammar::AbstractGrammar
)::Vector{Reaction}
    rule = grammar.rules[get_rule(program)]
    @match rule begin
        :(vcat(reaction,
            reaction_list)) => begin
            reaction::Reaction = interpret_reaction(program.children[1], grammar)
            reaction_list::Vector{Reaction} = interpret_network_reactions(
                program.children[2], grammar
            )
            return vcat(reaction, reaction_list)
        end

        :([reaction]) => begin
            reaction = interpret_reaction(program.children[1], grammar)
            return [reaction]
        end

        :(Vector{Reaction}()) => begin
            return Vector{Reaction}()
        end

        _ => throw(ArgumentError("Unknown rule: $rule"))
    end
end
