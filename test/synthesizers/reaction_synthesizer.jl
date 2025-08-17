@testitem "Reaction Synthesizer (from Atoms)" begin

    # Create some atoms
    a1 = Atom("H")
    a2 = Atom("O")

    # Synthesize options
    settings = SynthesizerSettings(max_depth = 8)
    reactions = synthesize_reactions([a1, a2], settings)
    @test length(reactions) > 0

    # TODO: Add more specific tests for the synthesized reactions
end
