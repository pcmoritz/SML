using DataStructures

# assumes that V[1:k] is sorted, and that element is in V[1:k]
function setdiff(V::Vector{Int}, element::Int, k::Int)
    index = 1
    while index <= k
        if V[index] == element
            break
        end
        index += 1
    end
    while index < k
        V[index] = V[index+1]
        index += 1
    end
end

function pendentpair!(func::Expr, keys::Vector{Float64}, used::Vector{Int}, V::Vector{Int}, k::Int)
    reset(func)
    fill!(used, 0)
    vnew = V[1]
    vold = -1
    used[vnew] = 1
    for i = 1:(k-1)
        vold = vnew
        incremental(func, vold)
        fill!(keys, Inf)
        for j = 1:k
            if used[j] == 1
                continue
            end
            keys[j] = incremental(func, V[j])
            reset(func, V[j])
        end
        argmin = -1
        curr_min = Inf
        for j = 1:n
            if keys[j] <= curr_min
                argmin = j
                curr_min = keys[j]
            end
        end
        vnew = V[argmin]
        used[argmin] = 1
    end
    return (vold, vnew)
end

function queyranne(func::Expr)
    p = size(func)
    s = zeros(p-1)
    G = DataStructures.IntDisjointSets(p)
    k = p
    keys = zeros(Float64, p)
    for i = 1:(p-1)
        (t, u) = pendentpair!(func, V, keys, k)
        # merge nodes u and t
        union!(G, t, u)
    end
end
