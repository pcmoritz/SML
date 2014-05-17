# (c) Philipp Moritz, 2014 and the Julia team

# This module contains a sorting method that doesn't allocate; see
# https://github.com/JuliaLang/SortingAlgorithms.jl

type RadixData
    vs::Vector{Int}
    ts::Vector{Int}
    bin::Array{Uint32,2}
    cbin::Vector{Uint32}
end

function RadixData(length::Int)
    vs = [1:length]
    ts = [1:length]
    bin = zeros(Uint32, 2^RADIX_SIZE, RADIX_ITERS)
    cbin = zeros(Uint32, 2^RADIX_SIZE)
    return RadixData(vs, ts, bin, cbin)
end

function mapping(x::Float64)
    y = reinterpret(Int64,x); uint64(y < 0 ? ~y : (y $ typemin(Int64)))
end

function cumsum!(from::Array{Uint32,2},to::Vector{Uint32},s::Uint32,i1::Int,n::Int,j::Int)
    if n < 128
        @inbounds to[i1] = s + from[i1,j]
        for i = i1+1:i1+n-1
            @inbounds to[i] = to[i-1] + from[i,j]
        end
    else
        n2 = div(n, 2)
        cumsum!(from, to, s, i1, n2,j)
        cumsum!(from, to, to[(i1+n2)-1], i1+n2, n-n2,j)
    end
end

const RADIX_SIZE = 11
const RADIX_MASK = 0x7FF
# Attention: This assumes that sizeof(Float64) == 8; if this is not
# the case, set RADIX_ITERS = iceil(sizeof(Float64)*8/RADIX_SIZE)
const RADIX_ITERS = 6

function radixsort!(X::Vector{Float64}, data::RadixData)
    lo = 1
    hi = length(X)

    # Init
    fill!(data.bin, 0)
    for i = 1:hi
        data.vs[i] = i
    end
    
    if lo > 1;  data.bin[1,:] = lo-1;  end

    # Histogram for each element, radix
    for i = lo:hi
        v = mapping(X[data.vs[i]])
        for j = 1:RADIX_ITERS
            idx = int((v >> (j-1)*RADIX_SIZE) & RADIX_MASK)+1
            @inbounds data.bin[idx,j] += 1
        end
    end

    # Sort!
    swaps = 0
    len = hi-lo+1
    for j = 1:RADIX_ITERS
        # Unroll first data iteration, check for degenerate case
        v = mapping(X[data.vs[hi]])
        idx = int((v >> (j-1)*RADIX_SIZE) & RADIX_MASK)+1

        # are all values the same at this radix?
        if data.bin[idx,j] == len;  continue;  end

        # cbin = cumsum(bin[:,j])
        cumsum!(data.bin, data.cbin, uint32(0), 1, 2^RADIX_SIZE, j)
        ci = data.cbin[idx]
        data.ts[ci] = data.vs[hi]
        data.cbin[idx] -= 1

        # Finish the loop...
        @inbounds for i in hi-1:-1:lo
            @inbounds v = mapping(X[data.vs[i]])
            idx = int((v >> (j-1)*RADIX_SIZE) & RADIX_MASK)+1
            @inbounds ci = data.cbin[idx]
            @inbounds data.ts[ci] = data.vs[i]
            @inbounds data.cbin[idx] -= 1
        end
        data.vs,data.ts = data.ts,data.vs
        swaps += 1
    end

    if isodd(swaps)
        data.vs,data.ts = data.ts,data.vs
        for i = lo:hi
            data.vs[i] = data.ts[i]
        end
    end
end
