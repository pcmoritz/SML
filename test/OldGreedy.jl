# (c) Philipp Moritz, 2014

using SortingAlgorithms

function polytope_corner(func, v)
    result = Array(Float64, length(v))
    SML.reset(func)
    min_val = typemax(Float64)
    acc = 0.0
    for i = 1:length(v)
        result[v[i]] = SML.incremental(func, v[i])
        acc += result[v[i]]
        min_val = min(min_val, acc)
    end
    return (result, min_val)
end

function greedy(func, w)
    v = sortperm(w)
    return polytope_corner(func, v)
end
