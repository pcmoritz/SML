# (c) Philipp Moritz, 2014

type IwataFunction <: Expr
    size :: Float64
    num :: Float64
end

IwataFunction(n) = IwataFunction(1.0 * n, 0)

function reset(func :: IwataFunction)
    func.num = 0.0
end

function incremental(func :: IwataFunction, element :: Int)
    val = - func.num * (func.size - func.num)
    val += (func.num + 1.0) * (func.size - func.num - 1.0)
    val -= (5.0 * element - 2.0 * func.size)
    func.num += 1.0
    return val
end

function evaluate(func :: IwataFunction, set :: Array{Int})
    num = sum(set)
    len = func.size
    result = 1.0 * num * (len - num)
    for i = 1:length(set)
      result -= (5.0 * i - 2.0 * len) * set[i]
    end
    return result;
end

emptyval(func::IwataFunction) = 0.0

variables(func :: IwataFunction) = begin
    return fill([1:iround(func.size)], 1)
end

curvature(func :: IwataFunction) = :submodular


