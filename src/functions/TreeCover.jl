# (c) Philipp Moritz, 2014

type TreeCoverFunction <: Expr
    d :: Int # height of the tree
    value :: Float64 # current value
    visited :: Vector{Bool}
end

TreeCoverFunction(d :: Int) = begin
    F = TreeCoverFunction(d, 0.0, fill(false, 2^d))
    reset(F)
    return F
end

reset(func :: TreeCoverFunction) = begin
    fill!(func.visited, false)
    func.value = 0.0
end

emptyval(func :: TreeCoverFunction) = 0.0

incremental(func :: TreeCoverFunction, element :: Int) = begin
    id = 1
    result = func.value
    # find most significant bit:
    number = element
    bitpos = 0
    while(number != 0)
        bitpos += 1
        number = number >> 1
    end
    # read each bit, beginning from the second most significant bit
    # and descent down the tree
    mask = 1 << (bitpos - 2)
    @inbounds for i = 1:bitpos-1
        if(!func.visited[id])
            func.visited[id] = true
            func.value += 1.0
        end
        bit = element & mask
        if bit == 0
            id = 2 * id
        else
            id = 2 * id + 1
        end
        element = element << 1
    end
    @inbounds if(!func.visited[id])
        func.visited[id] = true
        func.value += 1.0
    end
    return func.value - result
end
