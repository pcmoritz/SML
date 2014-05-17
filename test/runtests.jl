using SML
using Base.Test

reload("OldGreedy.jl")
reload("OldMinNormPoint.jl")

# Testing the Minimum Norm Point Algorithm (with Iwata's function)

function indicator(n, j)
    return [zeros(Int, j); ones(Int, n - j)]
end

f = SML.IwataFunction(1000)

@test SML.min_norm_point(f, randperm(1000)) == indicator(1000, 333)

# Testing cut functions

graph = SML.load_from_file("graph.dot")

cut_fn = SML.CutFunction(graph)

# Generated from A. Krause's toolbox with

# G=[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
# F_cut_fun = sfo_fn_cutfun(G);
# F_iw = sfo_fn_iwata(6);
# F = sfo_fn_lincomb([F_iw, F_cut_fun], [0.2, 1.0]);
# V = [1:6];
# A = sfo_min_norm_point(F,V) ==> [1 2 3 4 5 6]

F = 0.2 * SML.IwataFunction(6) + cut_fn

x = SML.evaluate(F, [1; 2; 3; 4; 5; 6])
y = SML.evaluate(F, [1; 4; 5])

A = SML.min_norm_point(F, [1:6], 1e-10)
a = SML.evaluate(F, SML.ind_to_set(A))

@test a <= x && a <= y


