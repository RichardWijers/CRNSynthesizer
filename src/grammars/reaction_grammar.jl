function reaction_grammar(;add_required_rule::Bool=false)
    grammar = @csgrammar begin
        reaction = Reaction(molecule_list, molecule_list)
        molecule_list = vcat([molecule], molecule_list)
        molecule_list = Vector{Molecule}()
    end

    if add_required_rule
        add_rule!(grammar, :(molecule_list = vcat([required_molecule], molecule_list)))
    end

    return grammar
end


function reaction_grammar(molecules::Vector{Molecule}; settings::SynthesizerSettings=SynthesizerSettings())
    grammar = reaction_grammar()

    for molecule in molecules
        add_rule!(grammar, :(molecule = $molecule))
    end

    # Make sure the reactions are atom balanced
    if !(haskey(settings.options, :disable_balanced_reaction) && settings.options[:disable_balanced_reaction])
        addconstraint!(grammar, BalancedReaction(complete_grammar=false))
    end

    # Makes sure that there is at least one molecule in the reaction
    addconstraint!(grammar, Forbidden(@c_rulenode 1{3, c}))
    addconstraint!(grammar, Forbidden(@c_rulenode 1{c, 3}))

    # Forbids reactions to be the same on both sides
    addconstraint!(grammar, Forbidden(@c_rulenode 1{a, a}))

    # Makes the molecule list ordered to break symmetries
    if !(haskey(settings.options, :disable_ordered_molecule_list) && settings.options[:disable_ordered_molecule_list])
        addconstraint!(grammar, Ordered((@c_rulenode 2{a, 2{b, c}}), [:a, :b]))
    end

    return grammar
end


function reaction_grammar(atoms::Vector{Atom}; settings::SynthesizerSettings=SynthesizerSettings())
    grammar = reaction_grammar()
    merge_grammars!(grammar, SMILES_grammar())

    for atom in atoms
        atom_str = "[" * string(atom) * "]"
        grammar = add_rule!(grammar, :(atom = $atom_str))
    end

    # Add the ValidSMILES constraint
    atom_dict, bond_dict = generate_atom_bond_dicts(grammar)
    digit_to_grammar = generate_digit_to_grammar(grammar)
    grammar_data = GrammarData(atom_dict, bond_dict, digit_to_grammar)
    addconstraint!(grammar, ValidSMILES(grammar_data))

    # Make sure the reactions are atom balanced
    if !(haskey(settings.options, :disable_balanced_reaction) && settings.options[:disable_balanced_reaction])
        addconstraint!(grammar, BalancedReaction(complete_grammar=true))
    end

    # Makes sure that there is at least one molecule in the reaction
    addconstraint!(grammar, Forbidden(@c_rulenode 1{3, c}))
    addconstraint!(grammar, Forbidden(@c_rulenode 1{c, 3}))

    # Forbids reactions to be the same on both sides
    addconstraint!(grammar, Forbidden(@c_rulenode 1{a, a}))

    # Makes the molecule list ordered to break symmetries
    if !(haskey(settings.options, :disable_ordered_molecule_list) && settings.options[:disable_ordered_molecule_list])
        addconstraint!(grammar, Ordered((@c_rulenode 2{a, 2{b, c}}), [:a, :b]))
    end
    
    return grammar
end