using CRNSynthesizer

using Test
using Aqua



@testset "Aqua" Aqua.test_all(CRNSynthesizer)

@testset "molecule.jl" include("datatypes/molecule.jl")
@testset "reaction_network.jl" include("datatypes/reaction_network.jl")