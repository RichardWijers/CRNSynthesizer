@testitem "Network Synthesizer (from Atoms)" begin

    # Create some atoms
    a1 = Atom("H")
    a2 = Atom("O")

    # Synthesize options
    settings = SynthesizerSettings(max_depth = 10)
    networks = synthesize_networks([a1, a2], settings)
    @test length(networks) > 0

    # TODO: Add more specific tests for the synthesized networks
end
