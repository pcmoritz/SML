# (c) Philipp Moritz, 2014

function polytope_corner!(func, v, result)
    reset(func)
    min_val = typemax(Float64)
    acc = 0.0
    for i = 1:length(v)
        result[v[i]] = incremental(func, v[i])
        acc += result[v[i]]
        min_val = min(min_val, acc)
    end
    return min_val
end

# ordering and result are passed by reference
function greedy(func, w, data::RadixData, result)
    radixsort!(w, data)
    polytope_corner!(func, data.vs, result)
end
