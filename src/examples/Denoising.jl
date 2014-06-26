D = 100
S = 30

img = zeros(Int, D, D)
img[S:(D-S), S:(D-S)] = 1

p = 0.2 # probability for pixel error
mask = rand(size(img)) .> 1-p

img[mask] = 1 .- img[mask]

ising = SML.IsingFunction(img, 1.0, 1.0, 1.0, 0.0)

A = SML.min_norm_point(ising, [1:D*D], 1e-12)
result = convert(Array{Int, 2}, reshape(A, (D, D)))

