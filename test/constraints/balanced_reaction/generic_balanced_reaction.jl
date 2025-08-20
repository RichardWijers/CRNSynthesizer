@testitem "GenericBalancedReaction (from Molecules)" begin
    using HerbConstraints
    using HerbCore

    # Create some molecules
    m1 = from_SMILES("[H]-[H]")
    m2 = from_SMILES("[O]=[O]")
    m3 = from_SMILES("[H]-[O]-[H]")

    # Create a reaction grammar
    grammar = reaction_grammar([m1, m2, m3])

    # Create a solver
    initial_node = @rulenode 1{2{4{5},2{6{7},3}},2{6{7},2{6{7},2{Hole[Bool[0, 0, 0, 0, 0, 1, 0, 1, 0]],Hole[Bool[0, 1, 1, 0, 0, 0, 0, 0, 0]]}}}}
    solver = GenericSolver(grammar, initial_node)
    solver.state.isfeasible = true

    constraint = first(solver.state.active_constraints)
    HerbConstraints.propagate!(solver, constraint)

    solver.state.isfeasible
end
