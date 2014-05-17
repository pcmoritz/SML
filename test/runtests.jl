using SML
using Base.Test

# Testing the Minimum Norm Point Algorithm (with Iwata's function)

function indicator(n, j)
    return [zeros(Int, j); ones(Int, n - j)]
end

f = SML.IwataFunction(1000)

@test SML.min_norm_point(f, randperm(1000)) == indicator(1000, 333)
