# function network_grammar()
#     return @csgrammar begin
#         network = ReactionNetwork(reaction_list)
#         reaction_list = vcat(reaction, reaction_list)
#         reaction_list = [reaction]
#     end
# end

function network_grammar()
    return @csgrammar begin
        network = ReactionNetwork(reaction_list)
        reaction_list = vcat(reaction, reaction_list)
        reaction_list = Vector{Reaction}()
    end
end



function network_grammar(atoms::Vector{Atom}; problem::ProblemDefinition=ProblemDefinition(), settings::SynthesizerSettings=SynthesizerSettings())

    grammar = network_grammar()
    merge_grammars!(grammar, reaction_grammar())
    merge_grammars!(grammar, SMILES_grammar())

    for atom in atoms
        atom_str = "[" * string(atom) * "]"
        add_rule!(grammar, :(atom = $atom_str))
    end

    # Add the ValidSMILES constraint
    atom_dict, bond_dict = generate_atom_bond_dicts(grammar)
    digit_to_grammar = generate_digit_to_grammar(grammar)
    grammar_data = GrammarData(atom_dict, bond_dict, digit_to_grammar)
    addconstraint!(grammar, ValidSMILES(grammar_data))

    # Make sure the reactions are atom balanced
    addconstraint!(grammar, BalancedReaction(complete_grammar=true))

    # Makes sure that there is at least one molecule in the reaction
    addconstraint!(grammar, Forbidden(@c_rulenode 4{6, c}))
    addconstraint!(grammar, Forbidden(@c_rulenode 4{c, 6}))

    # Forbids reactions to be the same on both sides
    addconstraint!(grammar, Forbidden(@c_rulenode 4{a, a}))

    # Makes the molecule list ordered to break symmetries
    addconstraint!(grammar, Ordered((@c_rulenode 5{a, 5{b, c}}), [:a, :b]))

    # Forbids reactions to occur more than once in the network
    addconstraint!(grammar, Forbidden(@c_rulenode 2{a, 2{a, c}}))

    # Makes the reaction list ordered to break symmetries
    if !(haskey(settings.options, :disable_ordered_reaction_list) && settings.options[:disable_ordered_reaction_list])
        addconstraint!(grammar, Ordered((@c_rulenode 2{a, 2{b, c}}), [:a, :b]))
    end

    # Make the ringbonds list ordered and forbid duplicates
    addconstraint!(grammar, Ordered((@c_rulenode 16{a, 16{b, c}}), [:a, :b]))
    addconstraint!(grammar, Forbidden((@c_rulenode 16{a, 16{a, c}})))


    required = problem.required_molecules
    required_rules = Vector{Tuple{Int, ReactionPosition}}()
    for (i, required_molecule) in enumerate(required)
        molecule = required_molecule.molecule
        add_rule!(grammar, :(required_molecule = $molecule))
        rule = findfirst(==(:($molecule)), grammar.rules)
        push!(required_rules, (rule, required_molecule.position))
    end


    if !(haskey(settings.options, :disable_contains_molecules) && settings.options[:disable_contains_molecules])
        add_rule!(grammar, :(molecule_list = vcat([required_molecule], molecule_list)))
        addconstraint!(grammar, ContainsMolecules(required_rules))
    end

    # println("Network grammar: ", grammar)

    return grammar
end

function network_grammar(molecules::Vector{Molecule}; problem::ProblemDefinition=ProblemDefinition(), required_molecules::Vector{Molecule}=Vector{Molecule}(), settings::SynthesizerSettings=SynthesizerSettings())
    grammar = network_grammar()
    merge_grammars!(grammar, reaction_grammar(add_required_rule=true))

    for molecule in molecules
        add_rule!(grammar, :(molecule = $molecule))
    end


    required = problem.required_molecules
    required_rules = Vector{Tuple{Int, ReactionPosition}}()
    for (i, required_molecule) in enumerate(required)
        molecule = required_molecule.molecule
        add_rule!(grammar, :(required_molecule = $molecule))
        rule = findfirst(==(:($molecule)), grammar.rules)
        push!(required_rules, (rule, required_molecule.position))
    end

    for molecule in required_molecules
        add_rule!(grammar, :(required_molecule = $molecule))
        rule = findfirst(==(:($molecule)), grammar.rules)
        push!(required_rules, (rule, UNKNOWN))
    end


    if length(required) == 0 && length(required_molecules) == 0
        addconstraint!(grammar, Forbidden(@c_rulenode 7))
    elseif !(haskey(settings.options, :disable_contains_molecules) && settings.options[:disable_contains_molecules])
        addconstraint!(grammar, ContainsMolecules(required_rules))
    end

    # Forbids reactions to occur more than once in the network
    addconstraint!(grammar, Forbidden(@c_rulenode 2{a, 2{a, c}}))

    if !(haskey(settings.options, :disable_ordered_reaction_list) && settings.options[:disable_ordered_reaction_list])
        # Makes the reaction list ordered to break symmetries
        addconstraint!(grammar, Ordered((@c_rulenode 2{a, 2{b, c}}), [:a, :b]))
    end

    # Forbids reactions to be the same on both sides
    addconstraint!(grammar, Forbidden(@c_rulenode 4{a, a}))


    # Make sure the reactions are atom balanced
    addconstraint!(grammar, BalancedReaction(complete_grammar=false, problem=problem))

    # # Makes sure that there is at least one molecule in the reaction
    addconstraint!(grammar, Forbidden(@c_rulenode 4{6, _}))
    addconstraint!(grammar, Forbidden(@c_rulenode 4{_, 6}))

    # Makes the molecule list ordered to break symmetries
    addconstraint!(grammar, Ordered((@c_rulenode 5{a, 5{b, _}}), [:a, :b]))
    addconstraint!(grammar, Ordered((@c_rulenode 7{a, 7{b, _}}), [:a, :b]))
    addconstraint!(grammar, Ordered((@c_rulenode 5{a, 7{b, _}}), [:a, :b]))
    addconstraint!(grammar, Ordered((@c_rulenode 7{a, 5{b, _}}), [:a, :b]))


    return grammar
end


function network_grammar(reactions::Vector{Reaction}; problem::ProblemDefinition=ProblemDefinition(), settings::SynthesizerSettings=SynthesizerSettings())
    grammar = network_grammar()

    for reaction in reactions[1:end-1]
        c_add_rule!(grammar, :(reaction = $reaction))
    end
    reaction = reactions[end]
    c_add_rule!(grammar, :(reaction = $reaction), recalculate=true)

    
    # Forbids reactions to occur more than once in the network
    addconstraint!(grammar, Forbidden(@c_rulenode 2{a, 2{a, b}}))

    if !(haskey(settings.options, :disable_ordered_reaction_list) && settings.options[:disable_ordered_reaction_list])
        # Makes the reaction list ordered to break symmetries
        addconstraint!(grammar, Ordered((@c_rulenode 2{a, 2{b, c}}), [:a, :b]))
    end


    reaction_to_rule = Dict{Any, Int}()
    for (i, rule) in enumerate(grammar.rules)
        reaction_to_rule[rule] = i
    end

    
    required_molecules::Dict{RequiredMolecule, Vector{Int}} = Dict{RequiredMolecule, Vector{Int}}()
    for required_molecule in problem.required_molecules
        required_molecules[required_molecule] = []

        for reaction in reactions
            if required_molecule.position == INPUT
                if any(input -> input[2] == required_molecule.molecule, reaction.inputs)
                    push!(required_molecules[required_molecule], reaction_to_rule[:($reaction)])
                end
            elseif required_molecule.position == OUTPUT
                if any(output -> output[2] == required_molecule.molecule, reaction.outputs)
                    push!(required_molecules[required_molecule], reaction_to_rule[:($reaction)])
                end
            else
                if any(input -> input[2] == required_molecule.molecule, reaction.inputs) ||
                   any(output -> output[2] == required_molecule.molecule, reaction.outputs)
                    push!(required_molecules[required_molecule], reaction_to_rule[:($reaction)])
                end
            end
        end
    end

    if !(haskey(settings.options, :disable_contains_molecules) && settings.options[:disable_contains_molecules])
        addconstraint!(grammar, ContainsReactions(required_molecules))
    end

    # println(grammar)

    return grammar
end



function network_grammar(problem::ProblemDefinition; settings::SynthesizerSettings=SynthesizerSettings())
    return network_grammar(get_atoms(problem); problem=problem, settings=settings)
end