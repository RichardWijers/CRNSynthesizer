@testitem "Molecule Synthesizer" begin

    # Create some atoms
    a1 = Atom("H")
    a2 = Atom("O")
    a3 = Atom("C")

    # Synthesize options
    settings = SynthesizerSettings(max_depth = 8)
    molecules = synthesize_molecules([a1, a2, a3], settings)
    @test length(molecules) > 0

    # TODO: Add more specific tests for the synthesized molecules
end
