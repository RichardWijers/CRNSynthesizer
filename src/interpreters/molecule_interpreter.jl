function interpret_molecule(program::AbstractRuleNode, grammar::AbstractGrammar)::Molecule
    rule = grammar.rules[get_rule(program)]

    if rule isa Molecule
        return rule
    end

    @match rule begin
        :(chain) => begin
            return from_SMILES(interpret_chain(program.children[1], grammar))
        end

        _ => begin
            child = program.children[1]
            mol = grammar.rules[get_rule(child)]
            if mol isa Molecule
                return mol
            else
                throw(ArgumentError("Unknown rule: $rule"))
            end
        end
    end
end

function interpret_chain(program::AbstractRuleNode, grammar::AbstractGrammar)::String
    rule = grammar.rules[get_rule(program)]

    @match rule begin
        :(SMILES_combine_chain(bond,
            structure,
            chain)) => begin
            bond_str = interpret_bond(program.children[1], grammar)
            structure_str = interpret_structure(program.children[2], grammar)
            chain_str = interpret_chain(program.children[3], grammar)
            return structure_str * bond_str * chain_str
        end

        :(structure * bond *
          chain) => begin
            structure_str = interpret_structure(program.children[1], grammar)
            bond_str = interpret_bond(program.children[2], grammar)
            chain_str = interpret_chain(program.children[3], grammar)
            return structure_str * bond_str * chain_str
        end

        :(atom *
          ringbonds) => begin
            atom_str = interpret_atom(program.children[1], grammar)
            ringbonds_str = interpret_ringbonds(program.children[2], grammar)
            return atom_str * ringbonds_str
        end

        _ => throw(ArgumentError("Unknown rule: $rule"))
    end
end

function interpret_bond(program::AbstractRuleNode, grammar::AbstractGrammar)::String
    return grammar.rules[get_rule(program)]
end

function interpret_structure(program::AbstractRuleNode, grammar::AbstractGrammar)::String
    rule = grammar.rules[get_rule(program)]

    @match rule begin
        :(atom * ringbonds *
          branches) => begin
            atom_str = interpret_atom(program.children[1], grammar)
            ringbonds_str = interpret_ringbonds(program.children[2], grammar)
            branches_str = interpret_branches(program.children[3], grammar)
            return atom_str * ringbonds_str * branches_str
        end

        _ => throw(ArgumentError("Unknown rule: $rule"))
    end
end

function interpret_atom(program::AbstractRuleNode, grammar::AbstractGrammar)::String
    return grammar.rules[get_rule(program)]
end

function interpret_ringbonds(program::AbstractRuleNode, grammar::AbstractGrammar)::String
    rule = grammar.rules[get_rule(program)]

    @match rule begin
        :(ringbond *
          ringbonds) => begin
            ringbond_str = interpret_ringbond(program.children[1], grammar)
            ringbonds_str = interpret_ringbonds(program.children[2], grammar)
            return ringbond_str * ringbonds_str
        end

        :("") => begin
            return ""
        end

        _ => throw(ArgumentError("Unknown rule: $rule"))
    end
end

function interpret_ringbond(program::AbstractRuleNode, grammar::AbstractGrammar)::String
    rule = grammar.rules[get_rule(program)]

    @match rule begin
        :(bond * digit) => begin
            bond_str = interpret_bond(program.children[1], grammar)
            digit_str = interpret_digit(program.children[2], grammar)
            return bond_str * digit_str
        end

        _ => throw(ArgumentError("Unknown rule: $rule"))
    end
end

function interpret_digit(program::AbstractRuleNode, grammar::AbstractGrammar)::String
    return grammar.rules[get_rule(program)]
end

function interpret_branches(program::AbstractRuleNode, grammar::AbstractGrammar)::String
    rule = grammar.rules[get_rule(program)]

    @match rule begin
        :(branch *
          branches) => begin
            branch_str = interpret_branch(program.children[1], grammar)
            branches_str = interpret_branches(program.children[2], grammar)
            return branch_str * branches_str
        end

        :("") => begin
            return ""
        end

        _ => throw(ArgumentError("Unknown rule: $rule"))
    end
end

function interpret_branch(program::AbstractRuleNode, grammar::AbstractGrammar)::String
    rule = grammar.rules[get_rule(program)]

    @match rule begin
        :("(" * bond * chain *
          ")") => begin
            bond_str = interpret_bond(program.children[1], grammar)
            chain_str = interpret_chain(program.children[2], grammar)
            return "(" * bond_str * chain_str * ")"
        end

        _ => throw(ArgumentError("Unknown rule: $rule"))
    end
end
