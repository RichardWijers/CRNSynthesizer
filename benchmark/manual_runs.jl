using CRNSynthesizer
using BenchmarkTools
using Term, Term.Progress
using Term: update!, with

include("data/estherification.jl")
include("data/water.jl")


pbar = ProgressBar()


# Simple problem
problem = water_problem()
settings = SynthesizerSettings(max_programs=1000, max_depth=10)

# Staged test
@time begin
    with(pbar) do
        global result = synthesize_networks(problem, settings, settings, initial_molecules_count=10, pbar=pbar)
    end
end

for network in result
    if !(length(network.reactions) == 1)
        continue
    end

    if !(length(get_molecules(network)) == 3)
        continue
    end

    println("Network: ", network)
end

# TODO: check if networks are equal


# Difficult problem
problem = estherification_problem(selected_known_indices=[2, 3, 5, 6], selected_expected_indices=[2, 3, 5, 6])
missing_information = estherification_problem(selected_known_indices=[1, 4], selected_expected_indices=[1, 4])

# Molecule test
settings = SynthesizerSettings(max_time=20, max_depth=8, benchmark_type=UntilFound, goal=get_molecules(missing_information))
atoms = get_atoms(problem)
molecules = synthesize_molecules(atoms, settings)


# Staged test
network = estherification_network()
molecule_settings = SynthesizerSettings(max_time=20, max_depth=8, benchmark_type=UntilFound, goal=get_molecules(missing_information))
network_settings = SynthesizerSettings(max_time=120, max_depth=8, benchmark_type=UntilFound, goal=Set([network]))
# @time begin
@profview begin
    with(pbar) do
        global result = synthesize_networks(problem, molecule_settings, network_settings, initial_molecules_count=30, pbar=pbar)
    end
end

# for network in result
# networks_subset = Set()
for network in result
    if !(length(network.reactions) == 2)
        continue
    end

    molecules = get_molecules(network)
    if !(length(molecules) == 6)
        continue
    end

    all_molecules = [
        from_SMILES("[H]-[C](=[O])-[O]-[H]"),
        from_SMILES("[C](-[C](-[H])(-[H])-[O]-[H])(-[H])(-[H])-[H]"),
        from_SMILES("[H]-[C](=[O])-[O]-[C](-[H])(-[H])-[C](-[H])(-[H])-[H]"),
        from_SMILES("[H]-[O]-[H]"),
        from_SMILES("[C](-[H])(-[H])(-[H])-[O]-[H]"),
        from_SMILES("[C](-[H])(=[O])-[O]-[C](-[H])(-[H])-[H]"),
    ]

    if !(issubset(all_molecules, molecules))
        continue
    end

    if compare_networks(network, estherification_network(), check_rates=false)
        println("Network: ", network)
        push!(networks_subset, network)
    end
end
