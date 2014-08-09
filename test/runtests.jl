# (c) Philipp Moritz, 2014

using SML
using Base.Test

cd(Pkg.dir("SML") * "/test/")

# Testing Iwata's function

iwata = SML.IwataFunction(100)
@test SML.currval(iwata) <= 1e-5
SML.incremental(iwata, 5)
@test (SML.currval(iwata) - 274.0) <= 1e-5
SML.incremental(iwata, 4)
@test (SML.currval(iwata) - 551.0) <= 1e-5
SML.reset(iwata, 5)
@test (SML.currval(iwata) - SML.evaluate(iwata, [4])) <= 1e-5

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

G = [0.0  1.0  1.2  0.0  0.0  0.0;
     1.3  0.0  1.4  0.0  0.0  0.0;
     1.5  1.6  0.0  1.7  0.0  0.0;
     0.0  0.0  1.8  0.0  1.9  2.0;
     0.0  0.0  0.0  2.1  0.0  2.2;
     0.0  0.0  0.0  2.3  2.4  0.0]

cut_fn = SML.CutFunction(SML.from_weight_matrix(G))

keys = zeros(Float64, SML.size(cut_fn))
used = zeros(Int, SML.size(cut_fn))

@test SML.curvature(cut_fn) == :submodular
@test SML.signature(cut_fn) == :pos

# More a test of the submodular ruleset:
@test SML.curvature(cut_fn * cut_fn) == :submodular

# Generated from A. Krause's toolbox with

# G=[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0;
#    0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
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

# Test composition of functions

p = 10

mu = abs(randn(p))
lambda = abs(randn(p))
T = 1000 * abs(randn(p))

lambda = lambda ./ (sum(lambda ./ mu) + 0.1)
rho = lambda ./ mu
z = rho .* T
y = rho ./ mu

h = SML.scalar_function(:increasing, :concave, x -> 1/(1.0 - x))
G = SML.Modular(z) + SML.Modular(-y) * SML.Composition(h, SML.Modular(rho))
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

A = SML.Modular([1, 2, 3])
B = SML.Modular([3, 4, 5])

F = A + B

@test SML.curvature(F) == :modular

C = SML.IwataFunction(3)

@test SML.curvature(A + C) == :submodular

# Test disciplined submodular programming

F = SML.Modular(rand(10))
@test SML.curvature(F * F) == :submodular

# Test Ising function

D = 100
S = 30

img = zeros(Int, D, D)
img[S:(D-S), S:(D-S)] = 1

p = 0.2 # probability for pixel error
mask = rand(size(img)) .> 1-p

img[mask] = 1 .- img[mask]

ising = SML.IsingFunction(img, 1.0, 1.0, 1.0, 0.0)

SML.evaluate_image(ising, img)

A = SML.min_norm_point(ising, [1:D*D], 1e-12)
result = convert(Array{Int, 2}, reshape(A, (D, D)))

@test abs(SML.emptyval(ising) + SML.incremental(ising, 1) + SML.incremental(ising, 2) + SML.incremental(ising, 12) - SML.evaluate(ising, [1, 2, 12])) <= 1e-4

# Testing GraphIntensity function

A = [0.0 1.0 1.0 0.0 0.0 0.0 0.0;
     0.0 0.0 1.0 1.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0 1.0 1.0 0.0;
     0.0 0.0 0.0 0.0 1.0 0.0 1.0;
     0.0 0.0 0.0 0.0 0.0 1.0 1.0;
     0.0 0.0 0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0 0.0 0.0 0.0]

g = SML.from_weight_matrix(A; is_directed=false)

intensity = SML.GraphIntensity(g)

@test abs(SML.incremental(intensity, 1) - 0.0) <= 1e-4
@test abs(SML.incremental(intensity, 2) - 1.0) <= 1e-4
@test abs(SML.incremental(intensity, 3) - 2.0) <= 1e-4
@test abs(SML.incremental(intensity, 4) - 1.0) <= 1e-4
@test abs(SML.incremental(intensity, 5) - 2.0) <= 1e-4

