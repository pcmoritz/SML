using SML
using Base.Test

cd(Pkg.dir("SML") * "/test/")

reload("OldGreedy.jl")
reload("OldMinNormPoint.jl")

# Testing the Minimum Norm Point Algorithm (with Iwata's function)

function indicator(n, j)
    return [zeros(Int, j); ones(Int, n - j)]
end

f = SML.IwataFunction(1000)

@test SML.min_norm_point(f, randperm(1000)) == indicator(1000, 333)

# Testing cut functions

# graph = SML.load_from_file("graph.dot")

G = [0.00000 1.00000 1.20000 0.00000 0.00000 0.00000;
     1.30000   0.00000   1.40000   0.00000   0.00000   0.00000;
     1.50000   1.60000   0.00000   1.70000   0.00000   0.00000;
     0.00000   0.00000   1.80000   0.00000   1.90000   2.00000;
     0.00000   0.00000   0.00000   2.10000   0.00000   2.20000;
     0.00000   0.00000   0.00000   2.30000   2.40000   0.00000]

cut_fn = SML.from_weight_matrix(G)

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

# Test Entropy Function

Sigma = [  2.26105    -0.71067    0.0398137  -0.448748;
 -0.71067     0.783202  -0.701331   -0.82428;
  0.0398137  -0.701331   0.937804    1.11777; 
 -0.448748   -0.82428    1.11777     2.10302]

F = SML.LogDet(Sigma)
@test abs(SML.incremental(F, 1) + SML.incremental(F, 2) - SML.evaluate(F, [1, 2])) <= 1e-4

# Test composition of functions

p = 10

mu = abs(randn(p))
lambda = abs(randn(p))
T = 1000 * abs(randn(p))

lambda = lambda ./ (sum(lambda ./ mu) + 0.1)
rho = lambda ./ mu
z = rho .* T
y = rho ./ mu

G = SML.Modular(z) + SML.Modular(-y) * SML.Composition(x -> 1/(1.0 - x), SML.Modular(rho))
A = SML.min_norm_point(G, [1:p], 1e-13)

a = SML.evaluate(G, SML.ind_to_set(A))

# Test exhaustively that we really found the minimum:

function powerset (x)
  result = {{}}
  for i in x, j = 1:length(result)
    push!(result, [result[j],i])
  end
  result
end

S = powerset(1:10)

for s in S
    X = Int[]
    for x in s
        X = [X; x]
    end
    b = SML.evaluate(G, X)
    @test a <= b
end


# Test set cover functions

F = SML.TreeCoverFunction(3)
G = SML.SetCoverFunction([Set(1), Set(1, 2), Set(1, 2, 4), Set(1, 2, 5), Set(1, 3), Set(1, 3, 6), Set(1, 3, 7)])

@test SML.evaluate(F, [1, 3]) == SML.evaluate(F, [1, 3])
