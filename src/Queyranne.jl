using DataStructures

# assumes that V[1:k] is sorted, and that element is in V[1:k]
function setdiff!(V::Vector{Int}, element::Int, k::Int)
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

function append!{T}(first::DLList{T}, second::DLList{T})
    first_len = length(first)
    second_len = length(second)
    second_last = second.prev
    first.prev.next = second.next
    second.prev.next = first.next
    second.next.prev = first.prev
    first.next.prev = second_last
    first.prev = second_last
    first.length = first_len + second_len
end

# first = DLList{Int}()
# second = DLList{Int}()
# push(first, 0)
# push(first, 1)
# push(second, 2)
# push(second, 3)
# append!(first, second)
# third = DLList{Int}()
# push(third, 4)
# append!(first, third)

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
    merged_fun = Induced(func, partitions, p, fill(0, p))
    keys = zeros(Float64, p)
    for i = 1:(p-1)
        (t, u) = pendentpair!(func, V, keys, k)
        # merge nodes u and t
        union!(G, t, u)
        setdiff!(V, u, k)
        k -= 1
    end
end
