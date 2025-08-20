@testitem "Reaction Synthesizer (from Atoms)" begin

    # Create some atoms
    a1 = Atom("H")
    a2 = Atom("O")

    # Synthesize options
    settings = SynthesizerSettings(max_depth = 8)
    reactions = synthesize_reactions([a1, a2], settings)
    @test length(reactions) > 0
    @test length(unique(reactions)) == length(reactions) broken=true
end


@testitem "Reaction Synthesizer (from Molecules)" begin

    # Create some molecules
    m1 = from_SMILES("[H]-[H]")
    m2 = from_SMILES("[O]=[O]")
    m3 = from_SMILES("[H]-[O]-[H]")

    # Synthesize options
    settings = SynthesizerSettings(max_depth = 5)
    reactions = synthesize_reactions([m1, m2, m3], settings)
    @test length(reactions) > 0
    @test length(unique(reactions)) == length(reactions)
end
