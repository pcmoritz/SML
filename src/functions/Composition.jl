# (c) Philipp Moritz, 2014

type Composition <: Expr
    phi::Function
    F::Expr
    accum::Float64
    last_val::Float64
end

Composition(phi::Function, F::Expr) = begin
    result = Composition(phi, F, 0.0, 0.0)
    reset(result)
    return result
end

emptyval(expr::Composition) = expr.phi(emptyval(expr.F))

reset(func::Composition) = begin
    reset(func.F)
    func.accum = emptyval(func.F)
    func.last_val = func.phi(func.accum)
end

incremental(func::Composition, element::Int) = begin
    r = incremental(func.F, element)
    func.accum += r
    temp = func.last_val
    func.last_val = func.phi(func.accum)
    return func.last_val - temp
end

evaluate(func::Composition, set) = begin
    return func.phi(evaluate(func.F, set))
end
