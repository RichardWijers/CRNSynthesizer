@enum BondType single double triple quadruple

function to_string(bond_type::BondType)::String
    return bond_type == single ? "-" :
           bond_type == double ? "=" : bond_type == triple ? "≡" : "≣"
end

struct Atom
    name::String
end

abstract type AbstractBond end

struct Bond <: AbstractBond
    from::Int
    to::Int
    bond_type::BondType

    function Bond(from::Int, to::Int, bond_type::BondType)
        from, to = minmax(from, to)
        return new(from, to, bond_type)
    end
end

struct Molecule
    atoms::Vector{Atom}
    bonds::Vector{AbstractBond}

    function Molecule(atoms::Vector{Atom}, bonds::Vector{Bond})
        # Create a copy of atoms and sort them
        sorted_indices = sortperm(atoms, by=atom -> atom.name)
        sorted_atoms = atoms[sorted_indices]
        
        # Create a mapping from old indices to new indices
        index_map = Dict(old_idx => new_idx for (new_idx, old_idx) in enumerate(sorted_indices))
        
        # Adjust bonds to reflect the new atom indices
        adjusted_bonds = map(bonds) do bond
            if bond isa Bond
                new_from = index_map[bond.from]
                new_to = index_map[bond.to]
                Bond(new_from, new_to, bond.bond_type)
            else
                bond  # For other types of bonds
            end
        end

        # Sort the bonds by the new atom indices
        sort!(adjusted_bonds, by=bond -> (bond.from, bond.to))
        
        return new(sorted_atoms, adjusted_bonds)
    end
end


# TODO: Think of a better Molecule struct that doesn't have this many symmetries
import Base.==
function ==(a::Molecule, b::Molecule)
    if length(a.atoms) != length(b.atoms) || length(a.bonds) != length(b.bonds)
        return false
    end

    if to_SMILES(a) != to_SMILES(b)
        return false
    end

    return true
end

import Base.hash
function hash(m::Molecule, h::UInt)
    h = hash(:Molecule, h)
    # println("Hashing molecule: ", m)
    h = hash(to_SMILES(m), h)
    return h
end

# import Base.hash
# function hash(m::Molecule, h::UInt)
#     h = hash(:Molecule, h)
#     # Hash sorted atom names
#     atom_names = Tuple(atom.name for atom in m.atoms)
#     h = hash(atom_names, h)
#     # Hash sorted bond tuples (from, to, bond_type as Int)
#     bond_tuples = Tuple((b.from, b.to, Int(b.bond_type)) for b in m.bonds if b isa Bond)
#     h = hash(bond_tuples, h)
#     return h
# end


Base.show(io::IO, atom::Atom) = print(io, atom.name)

function Base.show(io::IO, bond::Bond)
    print(io, "Bond(", bond.from, ", ", bond.to, ", ", to_string(bond.bond_type), ")")
end

Base.show(io::IO, molecule::Molecule) = begin
    print(io, "Molecule(atoms=[")
    print(io, join(molecule.atoms, ", "))
    print(io, "], bonds=[")
    print(io, join(molecule.bonds, ", "))
    print(io, "])")
end

function count_atoms(molecule::Molecule)::Dict{String, Int}
    atoms = Dict{String, Int}()
    for atom in molecule.atoms
        atoms[atom.name] = get(atoms, atom.name, 0) + 1
    end
    return atoms
end

function to_compact(molecule::Molecule)
    # Get the atom counts
    atoms = count_atoms(molecule)
    
    # Convert the compact representation to a string
    compact_str = ""
    for (atom_name, count) in atoms
        if count > 1
            compact_str *= "$atom_name$count"
        else
            compact_str *= atom_name
        end
    end

    return convert_to_subscript(compact_str)
end


function from_SMILES(smiles::String)
    # Extract all atoms
    atoms = Atom[]
    atom_matches = collect(eachmatch(r"\[(.*?)\]", smiles))
    for regmatch::RegexMatch in atom_matches
        atom_name = regmatch.captures[1]  # Get the content inside brackets
        push!(atoms, Atom(atom_name))
    end

    # Replace atoms with placeholders for easier parsing
    processed_smiles = smiles
    for (i, regmatch) in enumerate(atom_matches)
        submatch::String = Base.String(regmatch.match)
        processed_smiles = replace(processed_smiles, submatch => "[A$i]", count = 1)
    end

    current_atom_idx::Int = 0
    branch_stack = Int[]  # Stack to keep track of branching points
    ring_connections = Dict{Int, Tuple{Int, BondType}}()  # Store atom idx and bond type
    current_bond_type::BondType = single  # Default to single bond
    bonds = Bond[]  # Store bonds
    
    # Iterate safely through characters, handling Unicode properly
    i::Int = firstindex(processed_smiles)
    while i <= lastindex(processed_smiles)
        char = processed_smiles[i]

        if char == '['
            # Start of an atom placeholder
            atom_end = findnext(']', processed_smiles, i)
            if isnothing(atom_end)
                error("Unmatched '[' in SMILES string: $smiles")
            end

            atom_placeholder = processed_smiles[i:atom_end]
            atom_idx = parse(Int, match(r"A(\d+)", atom_placeholder).captures[1]::SubString{String})

            if current_atom_idx != 0  # If there's a previous atom, create a bond
                # Create bonds in both directions
                bond = Bond(current_atom_idx, atom_idx, current_bond_type)
                push!(bonds, bond)
            end

            current_atom_idx = atom_idx
            i = atom_end
        elseif char == '-'
            # Single bond (already default)
            current_bond_type = single
        elseif char == '='
            # Double bond
            current_bond_type = double
        elseif char == '#'
            # Triple bond
            current_bond_type = triple
        elseif char == '$'
            # Quadruple bond
            current_bond_type = quadruple
        elseif char == '≡'
            # Triple bond (Unicode)
            current_bond_type = triple
        elseif char == '≣'
            # Quadruple bond (Unicode)
            current_bond_type = quadruple
        elseif char == '('
            # Start of a branch
            push!(branch_stack, current_atom_idx)
        elseif char == ')'
            # End of a branch, return to parent atom
            if !isempty(branch_stack)
                current_atom_idx = pop!(branch_stack)
            end
        elseif isdigit(char)
            # Ring closure
            ring_num = parse(Int, string(char))

            if haskey(ring_connections, ring_num)
                # Close the ring
                other_atom_idx, ring_bond_type = ring_connections[ring_num]

                # Create bond
                bond = Bond(current_atom_idx, other_atom_idx, current_bond_type)
                push!(bonds, bond)

                # Reset bond type and remove the ring connection
                current_bond_type = single
                delete!(ring_connections, ring_num)
            else
                # Start the ring
                ring_connections[ring_num] = (current_atom_idx, current_bond_type)
                current_bond_type = single  # Reset for next bond
            end
        end
        i = nextind(processed_smiles, i)
    end

    return Molecule(atoms, bonds)
end








function to_SMILES(molecule::Molecule)::String
    if isempty(molecule.atoms)
        return ""
    end

    if length(molecule.atoms) == 1
        return "[" * molecule.atoms[1].name * "]"
    end

    # Create an adjacency and ringbond dict
    adjacency = Dict{Int, Vector{Tuple{Int, String}}}()
    ringbonds = Dict{Int, Vector{Tuple{Int, String}}}()
    ring_digit = 1

    visited_atoms = falses(length(molecule.atoms))
    visited_bonds = Dict{AbstractBond, Bool}()
    for bond in molecule.bonds
        visited_bonds[bond] = false
    end
    bonds_queue = Vector{Bond}()
    for bond in molecule.bonds
        if bond.from == 1
            push!(bonds_queue, bond)
        end
    end


    while !isempty(bonds_queue)
        bond = popfirst!(bonds_queue)
        visited_bonds[bond] = true

        # println("Processing bond: ", bond)

        if visited_atoms[bond.to] && visited_atoms[bond.from]
            if !haskey(ringbonds, bond.to)
                ringbonds[bond.to] = Vector{Int}()
            end
            if !haskey(ringbonds, bond.from)
                ringbonds[bond.from] = Vector{Int}()
            end
            push!(ringbonds[bond.to], (ring_digit, to_string(bond.bond_type)))
            push!(ringbonds[bond.from], (ring_digit, to_string(bond.bond_type)))
            ring_digit += 1
        else
            if !haskey(adjacency, bond.from)
                adjacency[bond.from] = Vector{Int}()
            end
            if !haskey(adjacency, bond.to)
                adjacency[bond.to] = Vector{Int}()
            end

            if visited_atoms[bond.to]
                push!(adjacency[bond.to], (bond.from, to_string(bond.bond_type)))
            else
                push!(adjacency[bond.from], (bond.to, to_string(bond.bond_type)))
            end

            visited_atoms[bond.from] = true
            visited_atoms[bond.to] = true

            for b in molecule.bonds
                if !(b in bonds_queue) && !visited_bonds[b] && (b.from == bond.to || b.to == bond.to || b.from == bond.from || b.to == bond.from)
                    push!(bonds_queue, b)
                end
            end
        end
    end

    function to_SMILES(atom_idx)
        result = "["* molecule.atoms[atom_idx].name *"]"

        if haskey(ringbonds, atom_idx)
            for ringbond in ringbonds[atom_idx]
                result *= ringbond[2] * string(ringbond[1])
            end
        end

        if !haskey(adjacency, atom_idx)
            return result
        end
        for (i, neighbour) in enumerate(adjacency[atom_idx])
            bond = neighbour[2]
            if i == length(adjacency[atom_idx])
                result *= bond* to_SMILES(neighbour[1])
            else
                result *= "(" *bond* to_SMILES(neighbour[1]) *")"
            end
        end

        return result
    end

    return to_SMILES(1)
end
















# function to_SMILES(molecule::Molecule)::String
#     if isempty(molecule.atoms)
#         return ""
#     end

#     # Create data structures for traversal
#     adjacency = Dict{Int, Vector{Tuple{Int, String}}}()
#     ringbonds = Dict{Int, Vector{Tuple{Int, String}}}()
#     visited = falses(length(molecule.atoms))
#     ring_digit = 1
    
#     # Build complete adjacency list from bonds
#     complete_adj = Dict{Int, Vector{Tuple{Int, BondType}}}()
#     for bond in molecule.bonds
#         if !haskey(complete_adj, bond.from)
#             complete_adj[bond.from] = []
#         end
#         if !haskey(complete_adj, bond.to)
#             complete_adj[bond.to] = []
#         end
#         push!(complete_adj[bond.from], (bond.to, bond.bond_type))
#         push!(complete_adj[bond.to], (bond.from, bond.bond_type))
#     end
    
#     # DFS function to build the traversal-based adjacency map
#     function dfs_build_adjacency(atom_idx, parent_idx=0)
#         visited[atom_idx] = true
        
#         if !haskey(adjacency, atom_idx)
#             adjacency[atom_idx] = []
#         end
        
#         # Process all neighbors
#         if haskey(complete_adj, atom_idx)
#             for (neighbor, bond_type) in complete_adj[atom_idx]
#                 # Skip parent we came from
#                 if neighbor == parent_idx
#                     continue
#                 end
                
#                 bond_str = to_string(bond_type)
                
#                 if visited[neighbor]
#                     # We found a ring closure
#                     if !haskey(ringbonds, atom_idx)
#                         ringbonds[atom_idx] = []
#                     end
#                     if !haskey(ringbonds, neighbor)
#                         ringbonds[neighbor] = []
#                     end
#                     push!(ringbonds[atom_idx], (ring_digit, bond_str))
#                     push!(ringbonds[neighbor], (ring_digit, bond_str))
#                     ring_digit += 1
#                 else
#                     # Normal bond, add to adjacency and continue DFS
#                     push!(adjacency[atom_idx], (neighbor, bond_str))
#                     dfs_build_adjacency(neighbor, atom_idx)
#                 end
#             end
#         end
#     end
    
#     # Start DFS from atom 1
#     dfs_build_adjacency(1)
    
#     # Recursive function to build the SMILES string (unchanged)
#     function to_SMILES(atom_idx)
#         result = "[" * molecule.atoms[atom_idx].name * "]"

#         if haskey(ringbonds, atom_idx)
#             for ringbond in ringbonds[atom_idx]
#                 result *= ringbond[2] * string(ringbond[1])
#             end
#         end

#         if !haskey(adjacency, atom_idx)
#             return result
#         end
        
#         for (i, neighbour) in enumerate(adjacency[atom_idx])
#             bond = neighbour[2]
#             if i == length(adjacency[atom_idx])
#                 result *= bond * to_SMILES(neighbour[1])
#             else
#                 result *= "(" * bond * to_SMILES(neighbour[1]) * ")"
#             end
#         end

#         return result
#     end

#     return to_SMILES(1)
# end