# (c) Philipp Moritz, 2014
# Test Cholesky updates and downdates and LogDet function

using Base.Test
import SML

n = 20

R = randn(n, n)
M = R' * R
S = Base.chol(M)
@test norm(S' * S - M) <= 1e-10

indices = BitArray(n)

k = 0
chol = zeros(n, n)

elements = [1, 2, 4]
s = length(elements)

for element = elements
    k = SML.chol_update!(M, indices, element, chol, k)
    indices[element] = true
end

@test norm(R[:,elements]' * R[:,elements] - (chol' * chol)[1:s,1:s]) <= 1e-10

elements = [1, 4]
s -= 1
SML.chol_downdate!(2, chol, 3)

@test norm(R[:,elements]' * R[:,elements] - (chol' * chol)[1:s,1:s]) <= 1e-10

# Test Entropy Function

Sigma = [  2.26105    -0.71067    0.0398137  -0.448748;
 -0.71067     0.783202  -0.701331   -0.82428;
  0.0398137  -0.701331   0.937804    1.11777; 
 -0.448748   -0.82428    1.11777     2.10302]

F = SML.LogDet(Sigma)
@test abs(SML.incremental(F, 1) + SML.incremental(F, 2) - SML.evaluate(F, [1, 2])) <= 1e-4
