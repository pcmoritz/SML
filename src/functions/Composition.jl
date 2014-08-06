# (c) Philipp Moritz, 2014

type Composition <: Expr
    phi::CVXExpr
    F # Vector{Expr}
    accum::Vector{Float64}
    last_val::Float64
end

# F is of type Vector{Expr} (can't encode that in Julia's type system)
Composition(phi::CVXExpr, F) = begin
    if size(phi) != length(F)
        error("Composition: Number of argument doesn't match arity of function")
    end
    outer_curvature = curvature(phi)
    for i = 1:size(phi)
        if outer_curvature[i] != :concave
            println("Disciplined Programming Error: Outer function is not concave in argument ", i)
        end
    end
    inner_curvature = map(curvature, F)
    inner_monotonicity = map(monotonicity, F)
    if composition(monotonicity(phi), inner_curvature, inner_monotonicity)
        result = Composition(phi, F, zeros(size(phi)), 0.0)
        reset(result)
        return result
    end
end

Composition(phi::CVXExpr, F::Expr) = begin
    return Composition(phi, [F])
end

emptyval(expr::Composition) = evaluate(expr.phi, map(emptyval, expr.F))

reset(func::Composition) = begin
    map(reset, func.F)
    func.accum = map(emptyval, func.F)
    func.last_val = evaluate(func.phi, func.accum)
end

incremental(func::Composition, element::Int) = begin
    n = size(func.phi)
    r = zeros(n)
    for i = 1:n
        r[i] = incremental(func.F[i], element)
    end
    func.accum += r
    temp = func.last_val
    func.last_val = evaluate(func.phi, func.accum)
    return func.last_val - temp
end

evaluate(func::Composition, set) = begin
    return evaluate(func.phi, map(x -> SML.evaluate(x, set), func.F))
end
