# (c) Philipp Moritz, 2014

type Modular <: Expr
    s::Array{Float64}
    Modular(x) = new(float64(x))
end

size(func::Modular) = Base.size(func.s, 1)

emptyval(func::Modular) = 0.0

evaluate(func::Modular, set) = begin
    result = 0.0
    for index = set
        result += func.s[index]
    end
    return result
end

incremental(func::Modular, element::Int) = begin
    return func.s[element]
end

reset(func :: Modular) = begin
end

variables(func :: Modular) = begin
    return singleton_partition(length(func.s))
end
