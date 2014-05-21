# (c) Philipp Moritz, 2014

using PyCall
using Graphs

typealias CutGraph Graph{ExVertex, ExEdge{ExVertex}}

type CutFunction <: Expr
    set :: Array{Int}
    val :: Float64
    graph :: CutGraph
end

CutFunction(g) = CutFunction(zeros(Int, num_vertices(g)), 0.0, g)

function from_weight_matrix(A)
    g = graph(ExVertex[], ExEdge{ExVertex}[]; is_directed=true)
    for x in 1:Base.size(A,1)
        add_vertex!(g, ExVertex(x, string(x)))
    end
    V = vertices(g)
    for j = 1:Base.size(A, 2)
        for i = 1:Base.size(A, 1)
            if A[i, j] != 0.0
                e = add_edge!(g, V[i], V[j])
                e.attributes["weight"] = A[i, j]
            end
        end
    end
    return CutFunction(g)
end

variables(func :: CutFunction) = begin
    G = func.graph
    V = map(vertex_index, vertices(G))
    E = map(e -> (vertex_index(source(e, G)),
                  vertex_index(target(e, G))), edges(G))
    g = simple_graph(length(V), is_directed=false)
    for (u, v) in E
        add_edge!(g, make_edge(g, u, v))
    end
    components = connected_components(g)
end

# reset the incremental oracle
function reset(func :: CutFunction)
    func.val = 0.0
    func.set = zeros(Int, length(func.set))
end

# evaluate one step of the incremental oracle
function incremental(func :: CutFunction, element :: Int)
    temp = func.val
    g = func.graph
    x = vertices(g)[element]
    for e in out_edges(x, g)
        v = target(e, g).index
        if func.set[v] == 0
            func.val += e.attributes["weight"]
        end
    end
    for e in in_edges(x, g)
        v = source(e, g).index
        if func.set[v] == 1
            func.val -= e.attributes["weight"]
        end
    end
    func.set[element] = 1
    return func.val - temp
end

# evaluate the whole function at once
function evaluate(func :: CutFunction, set :: Array{Int}; RESET=true)
    if RESET
        reset(func)
    end
    for elt in set
        incremental(func, elt)
    end
    return func.val
end

# If pydot is installed in your python distribution, you can uncomment
# the next line and use the following function:

# @pyimport pydot

# loading a graph from a graphviz file. we use the python library
# "pydot" for this purpose.
function load_from_file(filename)
    data = pydot.graph_from_dot_file(filename)
    cmd = "[(e.get_source(), e.get_destination(), e.get(\"weight\")) \
           for e in graph.get_edges()]"
    edges = pyeval(cmd, graph=data)
    cmd = "[v.get_name() for v in graph.get_nodes()]"
    nodes = pyeval(cmd, graph=data)
    g = graph(ExVertex[], ExEdge{ExVertex}[]; is_directed=true)
    for x in sort(map(int, nodes))
        add_vertex!(g, ExVertex(x+1, string(x+1)))
    end
    V = vertices(g)
    for (u, v, w) in edges
        e = add_edge!(g, V[int(u)+1], V[int(v)+1])
        e.attributes["weight"] = float(w)
    end
    return g
end
