# (c) Philipp Moritz, 2014
# Let G = (V, E, w) be a graph with weights; the cut function is given by
# F(X) = \sum_{i\in X, j\in V\X} w_{ij}

using PyCall
using Graphs

typealias CutGraph Graph{ExVertex, ExEdge{ExVertex}}

type CutFunction <: Expr
    set :: Array{Int}
    val :: Float64
    graph :: CutGraph
end

CutFunction(g) = CutFunction(zeros(Int, num_vertices(g)), 0.0, g)

size(func::CutFunction) = num_vertices(func.graph)

function from_weight_matrix(A; is_directed=true)
    g = graph(ExVertex[], ExEdge{ExVertex}[]; is_directed=is_directed)
    for x in 1:Base.size(A,1)
        add_vertex!(g, ExVertex(x, string(x)))
    end
    V = vertices(g)
    (I, J, val) = findnz(A)
    for i = 1:length(I)
        e = add_edge!(g, V[I[i]], V[J[i]])
        e.attributes["weight"] = val[i]
    end
    return g
end

function weight_signs(func::CutFunction)
    all_pos = true
    all_neg = true
    for e in Graphs.edges(func.graph)
        all_pos = all_pos && e.attributes["weight"] >= 0.0
        all_neg = all_neg && e.attributes["weight"] <= 0.0
    end
    return (all_pos, all_neg)
end

curvature(func::CutFunction) = begin
    (all_pos, all_neg) = weight_signs(func)
    if all_pos && all_neg
        return :modular
    end
    if all_pos
        return :submodular
    end
    if all_neg
        return :supmodular
    end
    return :unknown
end

monotonicity(func::CutFunction) = begin
    (all_pos, all_neg) = weight_signs(func)
    if all_pos
        return :increasing
    end
    if all_neg
        return :decreasing
    end
    return :unknown
end

signature(func::CutFunction) = begin
    (all_pos, all_neg) = weight_signs(func)
    if all_pos
        return :pos
    end
    if all_neg
        return :neg
    end
    return :unknown
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
function reset(func::CutFunction)
    func.val = 0.0
    fill!(func.set, 0)
end

emptyval(func::CutFunction) = 0.0
currval(func::CutFunction) = func.val

# evaluate one step of the incremental oracle
function incremental(func::CutFunction, element::Int)
    temp = func.val
    g = func.graph
    x = vertices(g)[element]
    func.set[element] = 1
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
    return func.val - temp
end

# undo the effect of incremental(func, element); only guaranteed to
# work if incremental(func, element) has been called before
function reset(func::CutFunction, element::Int)
    g = func.graph
    x = vertices(g)[element]
    for e in out_edges(x, g)
        v = target(e, g).index
        if func.set[v] == 0
            func.val -= e.attributes["weight"]
        end
    end
    for e in in_edges(x, g)
        v = source(e, g).index
        if func.set[v] == 1
            func.val += e.attributes["weight"]
        end
    end
    func.set[element] = 0
end

# evaluate the whole function at once
function evaluate(func::CutFunction, set::Array{Int}; RESET=true)
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
