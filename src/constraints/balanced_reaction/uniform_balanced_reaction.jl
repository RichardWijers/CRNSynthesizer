struct LocalUniformBalancedReaction <: AbstractLocalConstraint
    path::Vector{Int}
    input_atom_paths::Vector{Vector{Int}}
    output_atom_paths::Vector{Vector{Int}}
    input_molecule_paths::Vector{Vector{Int}}
    output_molecule_paths::Vector{Vector{Int}}
end

function get_atom_paths(solver::Solver, path::Vector{Int})
    node = get_node_at_location(solver, path)
    type = get_node_type(solver.grammar, node)
    if type == :atom
        return [path]
    end
    result = Vector{Vector{Int}}()
    for (i, child) in enumerate(node.children)
        next_path = push!(copy(path), i)
        result = vcat(result, get_atom_paths(solver, next_path))
    end
    return result
end

function get_molecule_paths(solver::Solver, path::Vector{Int})
    node = get_node_at_location(solver, path)
    type = get_node_type(solver.grammar, node)
    if type == :molecule || type == :required_molecule
        return [path]
    end

    result = Vector{Vector{Int}}()
    for (i, child) in enumerate(node.children)
        next_path = push!(copy(path), i)
        result = vcat(result, get_molecule_paths(solver, next_path))
    end
    return result
end

function LocalUniformBalancedReaction(solver, path::Vector{Int})
    # Find all input and output atom paths
    input_atom_paths = get_atom_paths(solver, child(path, 1))
    output_atom_paths = get_atom_paths(solver, child(path, 2))

    # Find all input and output molecule paths
    input_molecule_paths = get_molecule_paths(solver, child(path, 1))
    output_molecule_paths = get_molecule_paths(solver, child(path, 2))

    return LocalUniformBalancedReaction(
        path, input_atom_paths, output_atom_paths, input_molecule_paths, output_molecule_paths
    )
end

import Base.==
function ==(a::LocalUniformBalancedReaction, b::LocalUniformBalancedReaction)
    return a.path == b.path &&
           a.input_atom_paths == b.input_atom_paths &&
           a.output_atom_paths == b.output_atom_paths &&
           a.input_molecule_paths == b.input_molecule_paths &&
           a.output_molecule_paths == b.output_molecule_paths
end

import Base.hash
function hash(a::LocalUniformBalancedReaction, h::UInt)
    h = hash(a.path, h)
    h = hash(a.input_atom_paths, h)
    h = hash(a.output_atom_paths, h)
    h = hash(a.input_molecule_paths, h)
    h = hash(a.output_molecule_paths, h)
    return h
end

function HerbConstraints.shouldschedule(
        solver::Solver, constraint::LocalUniformBalancedReaction, path::Vector{Int}
)
    if path in constraint.input_atom_paths ||
       path in constraint.output_atom_paths ||
       path in constraint.input_molecule_paths ||
       path in constraint.output_molecule_paths
        return true
    end

    return false
end

function HerbConstraints.propagate!(
        solver::Solver, constraint::LocalUniformBalancedReaction
)
    # Helpers to parse and accumulate elemental counts
    # Extract element symbol from an atom rule value like "[H]"
    get_elem_from_atom_rule(rule_val) = begin
        s = String(rule_val)
        if startswith(s, "[") && endswith(s, "]")
            return s[2:end-1]
        else
            return s
        end
    end

    # Parse a count-encoded Symbol like :"C:1,H:2" into Dict{String,Int}
    function parse_count_symbol(x)::Dict{String, Int}
        d = Dict{String, Int}()
        s = String(x)
        if isempty(s)
            return d
        end
        parts = split(s, ',')
        for part in parts
            kv = split(part, ':')
            if length(kv) == 2
                d[kv[1]] = get(d, kv[1], 0) + parse(Int, kv[2])
            end
        end
        return d
    end

    function accumulate!(totals::Dict{String,Int}, add::Dict{String,Int})
        for (k,v) in add
            totals[k] = get(totals, k, 0) + v
        end
    end

    # Accumulate counts from atoms and molecules on one side
    function side_counts(atom_paths::Vector{Vector{Int}}, molecule_paths::Vector{Vector{Int}})
        counts = Dict{String, Int}()
        unfilled_atoms = 0
        unfilled_molecules = 0

        # Atom nodes (only count fixed ones)
        for p in atom_paths
            node = get_node_at_location(solver, p)
            if isfilled(node)
                rule_idx = HerbCore.get_rule(node)
                rule_val = solver.grammar.rules[rule_idx]
                elem = get_elem_from_atom_rule(rule_val)
                counts[elem] = get(counts, elem, 0) + 1
            else
                unfilled_atoms += 1
            end
        end

        # Molecule child nodes (count fixed molecules when determinable)
        for p in molecule_paths
            node = get_node_at_location(solver, p)
            rule_idx = HerbCore.get_rule(node)
            rule_val = solver.grammar.rules[rule_idx]
            if rule_val isa Molecule
                accumulate!(counts, count_atoms(rule_val))
            elseif rule_val isa Symbol
                accumulate!(counts, parse_count_symbol(rule_val))
            else
                println("Warning: Unexpected rule value type in molecule path $p: $rule_val")
                # Not directly countable (e.g., chain); we'll rely on atom nodes
                # present elsewhere, so do nothing here.
            end
        end

        return counts, unfilled_atoms
    end

    input_counts, unfilled_input_atoms = side_counts(
        constraint.input_atom_paths, constraint.input_molecule_paths
    )
    output_counts, unfilled_output_atoms = side_counts(
        constraint.output_atom_paths, constraint.output_molecule_paths
    )


    # println()
    # println("Input counts: ", input_counts)
    # println("Output counts: ", output_counts)
    # println("Unfilled input atoms: ", unfilled_input_atoms)
    # println("Unfilled output atoms: ", unfilled_output_atoms)



    # If both sides are fully filled, check if they match
    if unfilled_input_atoms == 0 && unfilled_output_atoms == 0 && !(input_counts == output_counts)
        HerbConstraints.set_infeasible!(solver)
        return
    end

    if unfilled_input_atoms == 0 && unfilled_output_atoms > 0
        # Check the counts difference
        difference = copy(input_counts)
        for (k, v) in output_counts
            difference[k] = v - get(input_counts, k, 0)
        end
        # println("Difference: ", difference)
    end
end
