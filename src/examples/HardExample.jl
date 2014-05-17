# (c) Philipp Moritz, 2014

# The example is from http://www.cs.berkeley.edu/~stefje/papers/smin.pdf

import SML

function hard(p::Int, M, W, G)
    q = ifloor(1.1 * p)

    # G = ifloor(rand(0:110, p, q)/100.0) # sparse biadjacency matrix

    G = G[1:p, 1:q]

    m = SML.Modular(M[1:p])
    w = SML.Modular(W[1:q])

    F = m + SML.Composition(x -> 100 * sqrt(x), SML.Induced(w, G))

    SML.min_norm_point(F, [1:p], 1e-12)
end

G = ifloor(rand(0:150, 10000, 10000)/100.0)

M = rand(-200:200, 10000)
W = rand(1:50, 10000)

history = [0 0.0]

for i = 1000:500:10000
    tic(); hard(i, M, W, G); time = toc();
    history = [history; i time]
end

for i = 1:size(history, 1)
    println("[", history[i,1], ", ", history[i,2], "],")
end
