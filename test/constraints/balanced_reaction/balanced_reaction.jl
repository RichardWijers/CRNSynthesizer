@testitem "BalancedReaction (from Molecules)" begin
    using HerbGrammar

    # Create some molecules
    m1 = from_SMILES("[H]-[H]")
    m2 = from_SMILES("[O]=[O]")
    m3 = from_SMILES("[H]-[O]-[H]")

    # Settings for the synthesizer without the BalancedReaction constraint
    settings = SynthesizerSettings(
        max_depth = 6, options = Dict{Symbol, Any}(:disable_balanced_reaction => true)
    )

    # Create a reaction grammar without the BalancedReaction constraint
    without_constraint_grammar = reaction_grammar([m1, m2, m3], settings = settings)

    # Create a reaction grammar with the BalancedReaction constraint
    with_constraint_grammar = reaction_grammar([m1, m2, m3], settings = settings)
    constraint = BalancedReaction()
    addconstraint!(with_constraint_grammar, constraint)
    with_constraint_grammar

    # Synthesize both options
    iterator = get_iterator(settings, without_constraint_grammar, :reaction)
    interpreter = x -> interpret_reaction(x, without_constraint_grammar)
    without_candidates = Vector{Reaction}()
    find_programs!(iterator, settings, interpreter, without_candidates)

    iterator = get_iterator(settings, with_constraint_grammar, :reaction)
    interpreter = x -> interpret_reaction(x, with_constraint_grammar)
    with_candidates = Vector{Reaction}()
    find_programs!(iterator, settings, interpreter, with_candidates)
    println()

    # Check the results
    @test length(without_candidates) > 0
    @test length(with_candidates) > 0
    @test length(with_candidates) < length(without_candidates)

    # Check the validity of the candidates
    valid_candidates = filter(x -> is_valid(x, constraint), without_candidates)
    @test length(valid_candidates) == length(with_candidates)
    @test all(x -> is_valid(x, constraint), with_candidates)
    @test all(x -> x in valid_candidates, with_candidates)
    @test all(x -> x in with_candidates, valid_candidates)
end
