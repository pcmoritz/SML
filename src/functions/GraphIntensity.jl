using Graphs

typealias SimpleGraph Graph{ExVertex, ExEdge{ExVertex}}

type GraphIntensity <: Expr
    set :: Array{Int}
    val :: Float64
    graph :: CutGraph
end

GraphIntensity(g) = GraphIntensity(zeros(Int, num_vertices(g)), 0.0, g)

function reset(func::GraphIntensity)
    func.val = 0.0
    fill!(func.set, 0)
end

emptyval(func::GraphIntensity) = 0.0

# assumes that element is not in func.set
function incremental(func::GraphIntensity, element::Int)
    temp = func.val
    g = func.graph
    x = vertices(g)[element]
    for e in out_edges(x, g)
        v = target(e, g).index
        if func.set[v] == 1
            func.val += e.attributes["weight"]
        end
    end
    func.set[element] = 1
    return func.val - temp
end
