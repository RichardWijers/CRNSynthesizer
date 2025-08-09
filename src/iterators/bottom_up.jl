@programiterator mutable BUIterator(bank) <: BottomUpIterator
function HerbSearch.init_combine_structure(iter::BUIterator)
    Dict(:max_combination_depth => 10)
end

function BUIterator(grammar::AbstractGrammar, start_symbol::Symbol; max_depth::Int = 10)
    return BUIterator(grammar, start_symbol, nothing; max_depth = max_depth)
end

function HerbSearch.combine(iter::BUIterator, state)
    addresses = Vector{CombineAddress}()
    max_in_bank = maximum(keys(HerbSearch.get_bank(iter)))
    max_per_type = Dict{Symbol, Int}()
    for (size, type_dict) in HerbSearch.get_bank(iter)
        for (type_sym, programs) in type_dict
            if !haskey(max_per_type, type_sym)
                max_per_type[type_sym] = 0
            end
            max_per_type[type_sym] = max(max_per_type[type_sym], size)
        end
    end

    grammar = get_grammar(iter.solver)
    terminals = grammar.isterminal
    nonterminals = .~terminals
    non_terminal_shapes = UniformHole.(partition(Hole(nonterminals), grammar), ([],))

    # if we have exceeded the maximum number of programs to generate
    if max_in_bank >= state[:max_combination_depth]
        return nothing, nothing
    end

    #check bound function
    function check_bound(type, combination)
        return 1 + sum((x[1] for x in combination)) > max_per_type[type]
    end

    function appropriately_typed(child_types)
        return combination -> child_types == [x[2] for x in combination]
    end

    # loop over groups of rules with the same arity and child types
    for shape in non_terminal_shapes
        type = HerbSearch.get_type(grammar, shape)
        child_types = grammar.childtypes[findfirst(shape.domain)]
        nchildren = length(child_types)

        # *Lazily* collect addresses, their combinations, and then filter them based on `check_bound`
        all_addresses = (
            (key, typename, idx) for key in keys(HerbSearch.get_bank(iter)) for
        typename in keys(HerbSearch.get_bank(iter)[key]) for
        idx in eachindex(HerbSearch.get_bank(iter)[key][typename])
        )
        all_combinations = Iterators.product(
            Iterators.repeated(all_addresses, nchildren)...
        )
        bounded_combinations = Iterators.filter(
            combination -> check_bound(type, combination), all_combinations
        )
        bounded_and_typed_combinations = Iterators.filter(
            appropriately_typed(child_types), bounded_combinations
        )

        # Construct the `CombineAddress`s from the filtered combinations
        append!(
            addresses,
            map(
                address_pair -> CombineAddress(shape, AccessAddress.(address_pair)),
                bounded_and_typed_combinations
            )
        )
    end

    return addresses, state
end
