@testitem "ReactionNetwork" begin
    # Test the ReactionNetwork constructor
    reaction1 = Reaction(1.0, [(2, from_SMILES("[H]-[H]"))], [(2, from_SMILES("[O]=[O]"))])
    reaction2 = Reaction(1.0, [(2, from_SMILES("[H]-[H]"))], [(2, from_SMILES("[O]=[O]"))])
    rn = ReactionNetwork([reaction1, reaction2])
    @test length(rn.reactions) == 2
end

@testitem "is_valid" begin
    # Test the is_valid function
    water = from_SMILES("[H]-[O]-[H]")
    hydrogen = from_SMILES("[H]-[H]")
    oxygen = from_SMILES("[O]=[O]")
    reaction = Reaction(1.0, [(2, hydrogen), (1, oxygen)], [(2, water)])
    rn = ReactionNetwork([reaction])
    @test is_valid(rn) == true
end
