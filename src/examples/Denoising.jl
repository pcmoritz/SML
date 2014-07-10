using Images
using ImageView

D = 100
S = 30

img = zeros(Int, D, D+10)
img[S:(D-S), S:(D+10-S)] = 1

p = 0.2 # probability for pixel error
mask = rand(size(img)) .> 1-p

img[mask] = 1 .- img[mask]

ImageView.display(Image(255*img))

ising = SML.IsingFunction(img, 1.0, 1.0, 1.0, 0.0)

A = SML.min_norm_point(ising, [1:D*(D+10)], 1e-12)
result = convert(Array{Int, 2}, reshape(A, (D, D+10)))

ImageView.display(Image(255*result))

# Dilbert example

cd(Pkg.dir("SML") * "/data/")

img = convert(Array, imread("denoise.ppm"))

# convert to 0/1
data = 1 * (img[:,:,3] .!= 0xff)
p = 0.1
mask = rand(size(data)) .> 1-p

new = copy(data)
new[mask] = 1 - new[mask]
new_img = Image(255*(1-new))
ImageView.display(new_img)

ising = SML.IsingFunction(new, 6.0, 1.0, 1.0, 1.0)

A = SML.min_norm_point(ising, [1:D*D], 1e-12)
result = convert(Array{Int, 2}, reshape(A, (D, D)))
img = Image(255*(1-result))
ImageView.display(img)
