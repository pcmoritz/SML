# (c) Philipp Moritz, 2014

abstract Expr

submodular(func) = begin
    curv = curvature(func)
    return (curv == :modular) || (curv == :submodular)
end

type Sum <: Expr
    first::Expr
    second::Expr
end

size(expr::Sum) = size(expr.first)

reset(expr::Sum) = begin
    reset(expr.first)
    reset(expr.second)
end

incremental(expr::Sum,element::Int) = 
    incremental(expr.first, element) + incremental(expr.second, element)

evaluate(expr::Sum, set) =
    evaluate(expr.first, set) + evaluate(expr.second, set)

curvature(expr::Sum) =
    plus(curvature(expr.first), curvature(expr.second))

+(first::Expr, second::Expr) = begin
    return Sum(first, second)
end

type Const <: Expr
    val::Float64
end

incremental(expr::Const, element::Int) = 0.0

evaluate(expr::Const, set) = expr.val

reset(expr::Const) = begin
end

emptyval(expr::Const) = expr.val

type Prod <: Expr
    first::Expr
    second::Expr
    accum::Float64
    firstval::Float64
    secondval::Float64
end

Prod(first::Expr, second::Expr) = begin
    result = Prod(first, second, 0.0, 0.0, 0.0)
    reset(result)
    return result
end

incremental(expr::Prod, element::Int) = begin
    first_d = incremental(expr.first, element)
    second_d = incremental(expr.second, element)
    incr = first_d*second_d + first_d*expr.secondval + second_d*expr.firstval
    expr.accum += incr
    expr.firstval += first_d
    expr.secondval += second_d
    return incr
end

evaluate(expr::Prod, set) = begin
    return evaluate(expr.first, set) * evaluate(expr.second, set)
end

reset(expr::Prod) = begin
    reset(expr.first)
    reset(expr.second)
    expr.firstval = emptyval(expr.first)
    expr.secondval = emptyval(expr.second)
    expr.accum = expr.firstval * expr.secondval
end

curvature(expr::Prod) = begin
    submod = true
    if !(submodular(expr.first) &&
         monotonicity(expr.first) == :increasing &&
         signature(expr.first) == :pos)
        submod = false
    end
    if !(submodular(expr.second) &&
         monotonicity(expr.second) == :increasing &&
         signature(expr.second) == :pos)
        submod = false
    end
    if submod
        return :submodular
    end
    return :unknown
end


*(first::Expr, second::Expr) = Prod(first, second)

*(first::Float64, second::Expr) = Prod(Const(first), second)

