using Test
using Aqua
using CRNSynthesizer



# @testset "Aqua" Aqua.test_all(
#     CRNSynthesizer;
#     stale_deps=(ignore=[:Aqua, :Test],)
# )

@testset "molecule.jl" include("datatypes/molecule.jl")
@testset "reaction_network.jl" include("datatypes/reaction_network.jl")