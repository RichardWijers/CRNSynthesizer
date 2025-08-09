function SMILES_combine_chain(bond, structure, chain)
    structure * bond * chain
end

function SMILES_grammar(
        atoms::Vector{Atom}; settings::SynthesizerSettings = SynthesizerSettings()
)
    grammar = @csgrammar begin
        molecule = chain

        chain = atom * ringbonds
        chain = SMILES_combine_chain(bond, structure, chain)
        # chain = structure * bond * chain

        structure = atom * ringbonds * branches

        branch = "(" * bond * chain * ")"
        branches = ""
        branches = branch * branches

        ringbond = bond * digit
        ringbonds = ""
        ringbonds = ringbond * ringbonds

        digit = "1" | "2" # | "3" | "4" | "5" | "6" | "7" | "8" | "9"

        bond = "-" | "=" | "≡" # | "≣"
    end

    for atom in atoms
        atom_str = "[" * string(atom) * "]"
        grammar = add_rule!(grammar, :(atom = $atom_str))
    end

    if !(
        haskey(settings.options, :disable_valid_smiles) &&
        settings.options[:disable_valid_smiles]
    )
        atom_dict, bond_dict = generate_atom_bond_dicts(grammar)
        digit_to_grammar = generate_digit_to_grammar(grammar)
        grammar_data = GrammarData(atom_dict, bond_dict, digit_to_grammar)
        addconstraint!(grammar, ValidSMILES(grammar_data))
    end

    # Make the ringbonds list tail ended
    addconstraint!(grammar, Ordered((@c_rulenode 10{a, 10{b, c}}), [:a, :b]))
    addconstraint!(grammar, Forbidden((@c_rulenode 10{a, 10{a, c}})))

    return grammar
end
