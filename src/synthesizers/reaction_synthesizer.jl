
function synthesize_reactions(
        atoms::Vector{Atom}, settings::SynthesizerSettings
)::Vector{Reaction}
    grammar = reaction_grammar(atoms; settings = settings)
    iterator = get_iterator(settings, grammar, :reaction)

    candidates = Vector{Reaction}()
    start_time = time()
    for program in iterator
        reaction = interpret_reaction(program, grammar)

        push!(candidates, reaction)

        if check_stop_condition(settings, start_time, candidates, reaction)
            break
        end
    end

    if haskey(settings.options, :unique_candidates) && settings.options[:unique_candidates]
        candidates = unique(candidates)
    end

    return candidates
end

function synthesize_reactions(
        molecules::Vector{Molecule}, settings::SynthesizerSettings
)::Vector{Reaction}
    grammar = reaction_grammar(molecules; settings = settings)
    iterator = get_iterator(settings, grammar, :reaction)

    # println(grammar)

    candidates = Vector{Reaction}()
    start_time = time()
    for program in iterator
        reaction = interpret_reaction(program, grammar)

        # println("\x1b[1mSynthesized reaction: $reaction\x1b[0m")

        push!(candidates, reaction)

        if check_stop_condition(settings, start_time, candidates, reaction)
            break
        end
    end

    if haskey(settings.options, :unique_candidates) && settings.options[:unique_candidates]
        candidates = unique(candidates)
    end

    return candidates
end
