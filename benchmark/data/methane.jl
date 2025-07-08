using CRNSynthesizer, Catalyst, OrdinaryDiffEq

function methane_problem(;
    selected_known_indices=1:4,
    selected_expected_indices=1:4)
    # Define the methane combustion reaction network
    rn = @reaction_network begin
        p1, CH₄ + 2O₂ --> CO₂ + 2H₂O
    end

    # Define the parameters
    tspan = (0.0, 10.0)
    u0 = [:CH₄ => 2.0, :O₂ => 4.0, :CO₂ => 0.0, :H₂O => 0.0]
    p = [:p1 => 0.2]

    # Solve the ODE problem
    prob = ODEProblem(rn, u0, tspan, p)
    sol = solve(prob, Tsit5())
    data_sol = solve(prob, Tsit5(), saveat=1.0)

    # Gather the time data and expected values
    time_data = data_sol.t[1:end]
    expected_CH4 = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:CH₄][1:end]
    expected_O2 = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:O₂][1:end]
    expected_CO2 = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:CO₂][1:end]
    expected_H2O = (0.9 .+ 0.2 * rand(length(time_data))) .* data_sol[:H₂O][1:end]

    # All possible molecules
    all_molecules = [
        from_SMILES("[H]-[C](-[H])(-[H])-[H]"),
        from_SMILES("[O]=[O]"),
        from_SMILES("[O]=[C]=[O]"),
        from_SMILES("[H]-[O]-[H]")
    ]
    
    # All expected profiles
    all_expected = [
        expected_CH4,
        expected_O2,
        expected_CO2,
        expected_H2O
    ]

    # Define the known molecules based on selected indices
    known_molecules = all_molecules[selected_known_indices]

    # Build expected profiles dictionary
    expected_profiles = Dict{Molecule, Vector{Float64}}()
    for (i, idx) in enumerate(selected_expected_indices)
        if idx in selected_known_indices
            known_idx = findfirst(j -> j == idx, selected_known_indices)
            expected_profiles[known_molecules[known_idx]] = all_expected[idx]
        end
    end

    # Define the problem
    problem = ProblemDefinition(
        known_molecules=known_molecules,
        expected_profiles=expected_profiles,
        time_data=time_data
    )

    return problem
end

function methane_network()
    # Define molecules using SMILES
    all_molecules = [
        from_SMILES("[H]-[C](-[H])(-[H])-[H]"),               # CH₄
        from_SMILES("[O]=[O]"),                               # O₂
        from_SMILES("[O]=[C]=[O]"),                           # CO₂
        from_SMILES("[H]-[O]-[H]")                            # H₂O
    ]

    # Define the reaction: CH₄ + 2O₂ --> CO₂ + 2H₂O
    reaction = CRNSynthesizer.Reaction(
        nothing,
        [(1, all_molecules[1]), (2, all_molecules[2])],
        [(1, all_molecules[3]), (2, all_molecules[4])]
    )

    return CRNSynthesizer.ReactionNetwork([reaction])
end
