# (c) Philipp Moritz, 2014

# Submodular function induced by a bipartite graph, see Stefanie
# Jegelka's PhD thesis "COMBINATORIAL PROBLEMS WITH SUBMODULAR
# COUPLING IN MACHINE LEARNING AND COMPUTER VISION", page 24

# The specified graph needs to be bipartite, with p nodes on the left;
# for the function we demand size(func) = q where q is the number of
# nodes on the right.

using Graphs

type Induced <: Expr
    func::Expr
    G::Vector{Vector{Int}}
    set::Vector{Int} # Nodes on the right side already accounted for
end

Induced(func::Expr, G::Array{Int, 2}) = begin
    p = Base.size(G, 1)
    q = Base.size(G, 2)
    graph = Array(Vector{Int}, p)
    for i = 1:p
        accum = Int[]
        for j = 1:q
            if G[i,j] == 1
                push!(accum, j)
            end
        end
        graph[i] = accum
    end
    return Induced(func, graph, fill(0, q))
end

size(expr::Induced) = Base.length(expr.G)

reset(expr::Induced) = begin
    fill!(expr.set, 0)
    reset(expr.func)
end

reset(expr::Induced, element::Int) = begin
    for next = reverse(expr.G[element])
        if expr.set[next] == 1
            reset(expr.func, next)
            expr.set[next] = 0
        end
    end
end

emptyval(expr::Induced) = emptyval(expr.func)
currval(expr::Induced) = currval(expr.func)

incremental(expr::Induced, element::Int) = begin
    acc = 0.0
    for next = expr.G[element]
        if expr.set[next] == 0
            # println("next ", next)
            # println("acc ", acc)
            # println("inner", expr.func)
            acc += incremental(expr.func, next)
            # println("acc ", acc)
            expr.set[next] = 1
        end
    end
    return acc
end
