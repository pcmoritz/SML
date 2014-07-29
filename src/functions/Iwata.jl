# (c) Philipp Moritz, 2014
# Iwata's function is defined via
# F(X) = |X||V\X| - \sum_{j\in X} (5j - 2n)

type IwataFunction <: Expr
    size::Float64
    num::Float64
    val::Float64
end

IwataFunction(n) = IwataFunction(1.0 * n, 0.0, 0.0)

function reset(func::IwataFunction)
    func.num = 0.0
end

function incremental(func::IwataFunction, element::Int)
    val = - func.num * (func.size - func.num)
    val += (func.num + 1.0) * (func.size - func.num - 1.0)
    val -= (5.0 * element - 2.0 * func.size)
    func.num += 1.0
    func.val += val
    return val
end

function reset(func::IwataFunction, element::Int)
    func.num -= 1.0
    val = - func.num * (func.size - func.num)
    val += (func.num + 1.0) * (func.size - func.num - 1.0)
    val -= (5.0 * element - 2.0 * func.size)
    func.val -= val    
    return nothing
end	 

function evaluate(func::IwataFunction, set::Array{Int})
    num = length(set)
    result = 1.0 * num * (func.size - num)
    for j = set
      result -= 5.0 * j - 2.0 * func.size
    end
    return result
end

emptyval(func::IwataFunction) = 0.0
currval(func::IwataFunction) = func.val

size(func::IwataFunction) = iround(func.size)

variables(func::IwataFunction) = begin
    return fill([1:iround(func.size)], 1)
end

curvature(func::IwataFunction) = :submodular

monotonicity(func::IwataFunction) = :unknown


