# (c) Philipp Moritz, 2014

type SetCoverFunction <: Expr
    sets :: Vector{Set{Int}}
    curr_set :: Vector{Bool}
    curr_val :: Float64
end

SetCoverFunction(sets::Vector{Set{Int}}) = begin
    p = maximum(map(maximum, sets))
    return SetCoverFunction(sets, fill(false, p), 0.0)
end

function reset(func::SetCoverFunction)
    func.curr_val = 0
    fill!(func.curr_set, false)
end

function incremental(func::SetCoverFunction, element::Int)
    result = 0
    for e in sets[element]
        if func.curr_set[e] == false
            result += 1
            func.curr_set[e] = true
        end
    end
    return result
end

function evaluate(func :: SetCoverFunction, set :: Array{Int}; RESET=true)
    if RESET
        reset(func)
    end
    for elt in set
        incremental(func, elt)
    end
    return func.curr_val
end

emptyval(func::SetCoverFunction) = 0.0
