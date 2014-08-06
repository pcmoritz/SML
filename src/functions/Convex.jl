# (c) Philipp Moritz, 2014
# A simple version of the convex ruleset

abstract CVXExpr

type ConvexFunction <: CVXExpr
    monotonicity::Vector{Symbol}
    curvature::Vector{Symbol}
    func::Function
    size::Int
end

square = ConvexFunction([:signed], [:convex], x -> x[1] * x[1], 1)
sqrt = ConvexFunction([:increasing], [:concave], x -> Base.sqrt(x[1]), 1)
log = ConvexFunction([:increasing], [:concave], x -> Base.log(x[1]), 1)
exp = ConvexFunction([:increasing], [:convex], x -> Base.log(x[1]), 1)
abs = ConvexFunction([:signed], [:convex], x -> Base.log(x[1]), 1)

scalar_function(mono::Symbol, curv::Symbol, func::Function) =
    ConvexFunction([mono], [curv], x -> func(x[1]), 1)

cobb_douglas(alpha::Vector{Float64}) = begin
    n = length(alpha)
    if all(alpha .>= 1.0)
        ConvexFunction([:decreasing for i in 1:n], [:concave for i in 1:n], x -> -prod(x .^ alpha), n)
    else
        error("Need to have alpha >= 1.")
    end
end

curvature(func::ConvexFunction) = func.curvature
monotonicity(func::ConvexFunction) = func.monotonicity
size(func::ConvexFunction) = func.size
evaluate(func::ConvexFunction, x::Vector{Float64}) = begin
    # if length(x) != 1
    #     error("Trying to evaluate a scalar function on a vector")
    # end
    # func.func(x[1])
    func.func(x)
end

type CVXSum <: CVXExpr
    first::CVXExpr
    second::CVXExpr
end

+(first::CVXExpr, second::CVXExpr) = begin
    if size(first) != size(second)
        throw(ArgumentError("+: Dimensions do not match"))
    end
    return CVXSum(first, second)
end

curvature(func::CVXSum) = begin
    result = Array(Symbol, size(func))
    first_curvature = curvature(func.first)
    second_curvature = curvature(func.second)
    for i = 1:size(func)
        result[i] = cvx_plus(first_curvature[i], second_curvature[i])
    end
    return result
end
evaluate(func::CVXSum, x::Vector{Float64}) = evaluate(func.first, x) + evaluate(func.second, x)
monotonicity(func::CVXSum) = begin
    result = Array(Symbol, size(func))
    first_monotonicity = monotonicity(func.first)
    second_monotonicity = monotonicity(func.second)
    for i = 1:size(func)
        result[i] = mon_plus(first_monotonicity[i], second_monotonicity[i])
    end
    return result
end
size(func::CVXSum) = size(func.first)

type CVXConst <: CVXExpr
    value::Float64
    size::Int
end

curvature(expr::CVXConst) = [:affine for i = 1:expr.size]
monotonicity(expr::CVXConst) = [:constant for i = 1:expr.size]
size(func::CVXConst) = func.size
evaluate(expr::CVXConst, x::Vector{Float64}) = expr.value

+(first::Real, second::CVXExpr) = begin
    return CVXConst(float(first), size(second)) + second
end

+(first::CVXExpr, second::Real) = second + first

type CVXProd <: CVXExpr
    first::Float64
    second::CVXExpr
end

evaluate(expr::CVXProd, x::Vector{Float64}) = expr.first * evaluate(expr.second, x)

curvature(expr::CVXProd) = begin
    curv = curvature(expr.second)
    result = Array(Symbol, size(expr))
    for i = 1:size(expr)
        result[i] = cvx_multiply(expr.first[i], curv[i])
    end
    return result
end

monotonicity(expr::CVXProd) = begin
    mon = monotonicity(expr.second)
    result = Array(Symbol, size(expr))
    for i = 1:size(expr.second)
        result[i] = mon_multiply(expr.first[i], mon[i])
    end
    return result
end
size(expr::CVXProd) = size(expr.second)

*(first::Real, second::CVXExpr) = CVXProd(float(first), second)
*(first::CVXExpr, second::Real) = second * first







