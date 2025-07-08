abstract type MyTopDown <: TopDownIterator end

@programiterator MyBFSIterator() <: MyTopDown

function priority_function(
    ::MyTopDown, 
    g::AbstractGrammar, 
    tree::AbstractRuleNode, 
    parent_value::Union{Real, Tuple{Vararg{Real}}},
    isrequeued::Bool
)
    #the default priority function is the bfs priority function
    if isrequeued
        return parent_value;
    end
    return parent_value + 1
end