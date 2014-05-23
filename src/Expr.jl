# (c) Philipp Moritz, 2014

abstract Expr

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


*(first::Expr, second::Expr) = Prod(first, second)

*(first::Float64, second::Expr) = Prod(Const(first), second)

