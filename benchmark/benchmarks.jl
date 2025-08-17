using BenchmarkTools
using CRNSynthesizer

include("data/estherification.jl")
include("data/water.jl")

const PROBLEMS = [
    (
        name = "Estherification",
        problem = estherification_problem(),
        max_iteration_goals = (molecules = 100, reactions = 100, networks = 100),
        until_found_goals = (
        # molecules = ,
        # reactions = ,
        # networks = ,
        )
    ),
    (
        name = "Water",
        problem = water_problem(),
        max_iteration_goals = (molecules = 100, reactions = 100, networks = 100),
        until_found_goals = (
        # molecules = ,
        # reactions = ,
        # networks = ,
        )
    )
]

const SUITE = BenchmarkGroup()

for p in PROBLEMS
    name = p.name
    problem = p.problem
    max_iteration_goals = p.max_iteration_goals
    until_found_goals = p.until_found_goals

    SUITE["$name - $(max_iteration_goals.molecules) molecules - from atoms"] = @benchmarkable begin
        atoms = get_atoms($problem)
        settings = SynthesizerSettings(
            max_programs = $max_iteration_goals.molecules, max_depth = 10
        )
        synthesize_molecules(atoms, settings)
    end

    # SUITE["$name - $(max_iteration_goals.reactions) reactions - from atoms"] = @benchmarkable begin
    #     atoms = get_atoms($problem)
    #     settings = SynthesizerSettings(max_programs=$max_iteration_goals.reactions, max_depth=10)
    #     synthesize_reactions(atoms, settings)
    # end

    # SUITE["$name - $(max_iteration_goals.reactions) reactions - from molecules"] = @benchmarkable begin
    #     molecules = get_molecules($problem)
    #     settings = SynthesizerSettings(max_programs=$max_iteration_goals.reactions, max_depth=10)
    #     synthesize_reactions(molecules, settings)
    # end

    # SUITE["$name - $(max_iteration_goals.networks) networks - from atoms"] = @benchmarkable begin
    #     atoms = get_atoms($problem)
    #     settings = SynthesizerSettings(max_programs=$max_iteration_goals.networks, max_depth=10)
    #     synthesize_networks(atoms, settings)
    # end

    # SUITE["$name - $(max_iteration_goals.networks) networks - from molecules"] = @benchmarkable begin
    #     molecules = synthesize_molecules(get_atoms($problem), SynthesizerSettings(max_programs=20, max_depth=10))
    #     settings = SynthesizerSettings(max_programs=$max_iteration_goals.networks, max_depth=10)
    #     synthesize_networks(molecules, settings)
    # end

    # SUITE["$name - $(max_iteration_goals.networks) networks - from reactions"] = @benchmarkable begin
    #     reactions = synthesize_reactions(get_atoms($problem), SynthesizerSettings(max_programs=100, max_depth=10))
    #     settings = SynthesizerSettings(max_programs=$max_iteration_goals.networks, max_depth=10)
    #     synthesize_networks(reactions, settings)
    # end

    # SUITE["$name - $(max_iteration_goals.networks) networks - problem -> networks"] = @benchmarkable begin
    #     settings = SynthesizerSettings(max_programs=$max_iteration_goals.networks, max_depth=10)
    #     synthesize_networks($problem, settings)
    # end

    # SUITE["$name - $(max_iteration_goals.networks) networks - problem -> molecules -> networks"] = @benchmarkable begin
    #     molecule_settings = SynthesizerSettings(max_depth=10)
    #     network_settings = SynthesizerSettings(max_programs=$max_iteration_goals.networks, max_depth=10)
    #     synthesize_networks($problem, molecule_settings, network_settings)
    # end

end

# Local usage:
# benchpkg CRNSynthesizer --rev=dirty --path=. --script=./benchmark/benchmarks.jl --add=HerbSearch#dev-bu-basic-interface
