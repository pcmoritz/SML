using DataStructures

function pendentpair!(func::Expr, keys::Vector{Float64}, singleton::Vector{Float64})
    reset(func)
    k = size(func)
    fill!(singleton, 0.0)
    vnew = 1
    vold = -1
    for i = 1:k
        singleton[i] = incremental(func, i)
        reset(func, i)
    end
    for i = 1:(k-1)
        vold = vnew
        prev_val = incremental(func, i)
        fill!(keys, Inf)
        for j = i:k
            incremental(func, j)
            keys[j] = currval(func) - singleton[j]
            print("keys ")
            print(j)
            print(" ")
            print(keys[j])
            print("\n")
            reset(func, j)
        end
        argmin = -1
        curr_min = Inf
        for j = i:k
            if keys[j] <= curr_min
                argmin = j
                curr_min = keys[j]
            end
        end
        vnew = argmin
    end
    return (vold, vnew)
end

function merge!(partitions, j, k)
    n = Base.length(partitions)
    partitions[j] = [partitions[j]; partitions[k]]
    partitions[k] = partitions[n]
    return partitions[1:n-1]
end

# partitions = merge!(partitions, 1, 20)
# partitions = merge!(partitions, 2, 3)
# partitions = merge!(partitions, 1, 2)
# partitions = merge!(partitions, 2, 3)
# partitions = merge!(partitions, 3, 4)
# partitions = merge!(partitions, 2, 3)

function queyranne(func::Expr)
    p = size(func)
    k = p
    partitions = Array(Vector{Int}, p)
    for i = 1:p
        partitions[i] = [i]
    end
    merged_fun = Induced(func, partitions, fill(0, p))
    singleton = zeros(Float64, p)
    keys = zeros(Float64, p)
    min_val = Inf
    min_set = [-1]
    for i = 1:(p-1)
        (t, u) = pendentpair!(merged_fun, keys, singleton)
        println("t ", t, " u ", u)
        reset(merged_fun)
        val = incremental(merged_fun, u)
        reset(merged_fun, u)
        if val <= min_val
            min_val = val
            min_set = partitions[u]
        end
        partitions = merge!(partitions, t, u)
        merged_fun.G = partitions
    end
    println(min_val)
    println(min_set)
end
