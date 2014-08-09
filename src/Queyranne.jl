using DataStructures

function pendentpair!(func::Expr, keys::Vector{Float64}, singleton::Vector{Float64}, used::BitArray)
    reset(func)
    k = size(func)
    fill!(singleton, 0.0)
    fill!(used, false)
    used[1] = true
    vnew = 1
    vold = -1
    for i = 1:k
        # println("i: ", i)
        singleton[i] = incremental(func, i)
        # println("func.set is ", func.set)
        # println("func.G is ", func.G)
        # println("value is ", singleton[i])
        reset(func, i)
    end
    # println("singletons ", singleton)
    for i = 1:(k-1)
        # println("i ", i)
        # println("vnew ", vnew)
        # println("vold ", vold)
        vold = vnew
        prev_val = incremental(func, vold)
        # println("prev_function", func.set)
        fill!(keys, Inf)
        for j = 1:k
            if used[j]
                continue
            end
            incremental(func, j)
            # println("partitions ", func.G)
            # println("post_function", func.set)
            keys[j] = currval(func) - singleton[j]
            # println("keys[j] ", j, " ", keys[j])
            reset(func, j)
        end
        argmin = -1
        curr_min = Inf
        for j = 1:k
            if keys[j] <= curr_min
                argmin = j
                curr_min = keys[j]
            end
        end
        used[argmin] = true
        vnew = argmin
        # println("vnew = argmin = ", argmin)
    end
    return (vold, vnew)
end

function merge!(partitions, j, k)
    n = Base.length(partitions)
    partitions[j] = [partitions[j]; partitions[k]]
    partitions[k] = partitions[n]
    return partitions[1:n-1]
end

function queyranne(func::Expr)
    p = size(func)
    k = p
    partitions = Array(Vector{Int}, p)
    for i = 1:p
        partitions[i] = [i]
    end
    merged_fun = Induced(func, partitions, fill(0, p))
    singleton = zeros(Float64, p)
    used = BitArray(p)
    keys = zeros(Float64, p)
    A = Array(Vector{Int}, p)
    z = zeros(p)
    for i = 1:p-1
        (t, u) = pendentpair!(merged_fun, keys, singleton, used)
        # println("pair ", t, u)
        reset(merged_fun)
        z[i] = incremental(merged_fun, u)
        A[i] = partitions[u]
        reset(merged_fun, u)
        partitions = merge!(partitions, t, u)
        # println("partitions ", partitions)
        merged_fun.G = partitions
    end
    min = Inf
    argmin = -1
    for i = 1:p-1
        if z[i] <= min
            min = z[i]
            argmin = i
        end
    end
    return A[argmin], min
end
