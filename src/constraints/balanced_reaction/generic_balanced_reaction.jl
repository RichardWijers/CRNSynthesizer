struct LocalGenericBalancedReaction <: AbstractLocalConstraint 
    path::Vector{Int}  # Path to the reaction node
end

import Base.==
function ==(a::LocalGenericBalancedReaction, b::LocalGenericBalancedReaction)
    return a.path == b.path
end

import Base.hash
function hash(a::LocalGenericBalancedReaction, h::UInt)
    return hash(a.path, h)
end


function HerbConstraints.shouldschedule(solver::Solver, constraint::LocalGenericBalancedReaction, path::Vector{Int})
    # Check if the path is a child of the reaction node
    if length(path) <= length(constraint.path)
        return false
    end
    
    if path[1:length(constraint.path)] != constraint.path
        return false
    end
    
    # Schedule for atom nodes and relevant holes
    node = get_node_at_location(solver, path)
    type = get_node_type(solver.grammar, node)

    # Only schedule for atoms or holes that could contain atoms
    return type == :atom || (node isa Hole && 
        (type == :chain || type == :structure || type == :branch || 
         type == :branches || type == :molecule_list || type == :molecule))
end


function get_relevant_paths(solver::Solver, c::LocalGenericBalancedReaction, path::Vector{Int})
    # Get the relevant paths for the current node
    # println(get_tree(solver))
    # println("path = $path")
    node = get_node_at_location(solver, path)
    type = get_node_type(solver.grammar, node)

    if isuniform(node)
        @match type begin
            :molecule_list => begin

                rule = HerbConstraints.get_rule(node)
                @match solver.grammar.rules[rule] begin
                    :(vcat([molecule], molecule_list)) => begin
                        # Get the relevant paths for the molecule and list
                        molecule_paths, molecule_holes, required_molecule_paths = get_relevant_paths(solver, c, push!(copy(path), 1))
                        list_paths, list_holes, required_molecule_list_paths = get_relevant_paths(solver, c, push!(copy(path), 2))

                        # Combine the paths and holes
                        return vcat(molecule_paths, list_paths), vcat(molecule_holes, list_holes), vcat(required_molecule_paths, required_molecule_list_paths)
                    end
                    :(vcat([required_molecule], molecule_list)) => begin
                        list_paths, list_holes, required_molecule_list_paths = get_relevant_paths(solver, c, push!(copy(path), 2))

                        return list_paths, list_holes, vcat(required_molecule_list_paths, push!(copy(path), 1))
                    end
                    :(Vector{Molecule}()) => return [], [], []
                    x => throw("Unknown molecule_list rule: $x")
                end
            end
            :molecule => begin
                return get_relevant_paths(solver, c, push!(copy(path), 1))
            end
            :chain => begin
                rule = HerbCore.get_rule(node)
                @match solver.grammar.rules[rule] begin
                    :(SMILES_combine_chain(bond, structure, chain)) => begin
                        # Get the relevant paths for the bond and structure
                        structure_paths, structure_holes, required_structure_paths = get_relevant_paths(solver, c, push!(copy(path), 2))
                        chain_paths, chain_holes, required_chain_paths = get_relevant_paths(solver, c, push!(copy(path), 3))

                        # Combine the paths and holes
                        return vcat(chain_paths, structure_paths), vcat(chain_holes, structure_holes), vcat(required_chain_paths, required_structure_paths)

                    end
                    :(atom * ringbonds) => begin
                        atom_path = push!(copy(path), 1)
                        return [atom_path], [], []
                    end
                    x => throw("Unknown chain rule: $x")
                end
            end
            :structure => begin
                atom_path = push!(copy(path), 1)
                branches_paths, branches_holes, required_branches_paths = get_relevant_paths(solver, c, push!(copy(path), 3))
                # Combine the paths and holes
                return vcat([atom_path], branches_paths), branches_holes, required_branches_paths
            end
            :ringbonds => begin
                return [], [], []
            end
            :branches => begin
                rule = HerbCore.get_rule(node)
                @match solver.grammar.rules[rule] begin
                    :(branch * branches) => begin
                        # Get the relevant paths for the branch and branches
                        branch_paths, branch_holes, required_branch_paths = get_relevant_paths(solver, c, push!(copy(path), 1))
                        branches_paths, branches_holes, required_branches_paths = get_relevant_paths(solver, c, push!(copy(path), 2))

                        # Combine the paths and holes
                        return vcat(branch_paths, branches_paths), vcat(branch_holes, branches_holes), vcat(required_branch_paths, required_branches_paths)
                    end
                    :("") => return [], [], []
                    x => throw("Unknown branches rule: $x")
                end
            end
            :branch => begin
                return get_relevant_paths(solver, c, push!(copy(path), 2))
            end
            x => throw("Unknown type: $x")
        end
    else
        @match type begin
            :molecule_list => begin
                return Vector{Vector{Int}}(), [path], []
            end
            :ringbonds => begin
                return [], [], []
            end
            :branches => begin
                return [], [path], []
            end
            :chain => begin
                return [], [path], []
            end
            :molecule => begin
                return [], [path], []
            end
            x => throw("Unknown type: $x")
        end
    end
end


function HerbConstraints.propagate!(solver::Solver, constraint::LocalGenericBalancedReaction)  
    node = get_node_at_location(solver, constraint.path)
    type = get_node_type(solver.grammar, node)
    
    if !isuniform(node)
        return
    end

    left_path = push!(copy(constraint.path), 1)
    right_path = push!(copy(constraint.path), 2)

    left_atoms, left_holes, required_left_molecules = get_relevant_paths(solver, constraint, left_path)
    right_atoms, right_holes, required_right_molecules = get_relevant_paths(solver, constraint, right_path)

    # println("Counting atoms and holes...")
    # println("Left side:")
    # println("Atoms: $(length(left_atoms)), Holes: $(length(left_holes))")
    # println("Right side:")
    # println("Atoms: $(length(right_atoms)), Holes: $(length(right_holes))")

    hole_terminals = findall(solver.grammar.isterminal)


    # Fix all the holes on the right side if there can't be more atoms
    if length(left_holes) == 0 && length(left_atoms) == length(right_atoms) && length(right_holes) > 0
        
        for hole in right_holes
            node = get_node_at_location(solver, hole)
            # println("Fixing hole: $node")
            c_remove_all_but!(solver, hole, hole_terminals)
            node = get_node_at_location(solver, hole)
            # println("Fixed hole: $node")
        end
        
    end

    # Both sides can not get any more atoms, so they should be the same
    if length(left_holes) == 0 && length(right_holes) == 0
        if length(left_atoms) != length(right_atoms)
            HerbConstraints.set_infeasible!(solver)
            # println("infeasible")
            return
        end
    end

    # Only the right side can get atoms, so it should be the same or more as the left side
    if length(left_holes) == 0 && length(right_holes) > 0
        if length(left_atoms) < length(right_atoms)
            HerbConstraints.set_infeasible!(solver)
            # println("infeasible")
            return
        end
    end

    # Only the left side can get atoms, so it should be the same or more as the right side
    if length(right_holes) == 0 && length(left_holes) > 0
        if length(right_atoms) < length(left_atoms)
            HerbConstraints.set_infeasible!(solver)
            # println("infeasible")
            return
        end
    end

    if length(left_holes) == 0 && length(right_holes) == 0
        # println("deactivated")
        HerbConstraints.deactivate!(solver, constraint)
        return
    end
end
