type CutFunctionTest <: Expr
    set :: Array{Int}
    G :: Array{Float64, 2}
end

CutFunctionTest(G) = CutFunctionTest(Array(Int, 0), G)

size(func::CutFunctionTest) = size(func.G, 1)

function incremental(func::CutFunctionTest, element::Int)
    new_set = [func.set, element]
    inc = SML.evaluate(func, new_set) - SML.evaluate(func, func.set)
    func.set = new_set
    return inc
end

function reset(func::CutFunctionTest, element::Int)
    func.set = setdiff(func.set, [element])
    return nothing
end

function reset(func::CutFunctionTest)
    func.set = Array(Int, 0)
    return nothing
end

function evaluate(func::CutFunctionTest, set::Array{Int})
    n = size(func)
    set_c = setdiff([1:n], set)
    return sum(func.G[set, set_c])
end

emptyval(func::CutFunctionTest) = 0.0
currval(func::CutFunctionTest) = SML.evaluate(func, func.set)

