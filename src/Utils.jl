using DataStructures

# Union of two sets represented by bitsets
function set_union(A :: Array{Int}, B :: Array{Int})
    result = zeros(Int, max(length(A), length(B)))
    for i = 1:length(result)
        if A[i] == 1 || B[i] == 1
            result[i] = 1
        end
    end
    return result
end

function create_partition_dict(partition)
    dict = Dict{Int, Int}()
    for i = 1:length(partition)
        for element = partition[i]
            dict[element] = i
        end
    end
    return dict
end

function merge_partitions!(partition :: Dict{Int, Int}, a :: Int, b :: Int)
    new_index = min(a, b)
    for (k, v) in partition
        if v == a || v == b
            partition[k] = new_index
        end
    end
    return partition
end

# Given two partitions of the numbers 1, \dots, n, find the finest
# possible partition that is coarser than both of them. Can likely be
# optimized from O(n^2) to O(n) if needed.
function align_partitions(a, b)
    dict_a = create_partition_dict(a)
    for set in b
        for elt in set
            merge_partitions!(dict_a, dict_a[set[1]], dict_a[elt])
        end
    end
    size = maximum(values(dict_a))
    result = Array(Array{Int}, size)
    for i = 1:size
        result[i] = Array(Int, 0)
    end
    for (k, v) in dict_a
        push!(result[v], k)
    end
    return result
end

# The partition of the set [1:n] into all singletons
function singleton_partition(n)
    result = Array(Array{Int}, n)
    for i = 1:n
        result[i] = [i]
    end
    return result
end

# Convert an indicator function 1_A to the set A
function ind_to_set(ind::Vector{Float64})
    A = Int[]
    for i = 1:length(ind)
        if ind[i] > 0.5
            push!(A, i)
        end
    end
    return A
end
