using Images
using ImageView
using SML

data_dir = ENV["DATA_DIR"]
# directory = joinpath(data_dir, "small-features")
# files = readdir(directory)

# data = convert(Array, imread(joinpath(directory, files[4])))

# data = convert(Array, imread(joinpath(data_dir, "features/features-small.png")))

data = convert(Array, imread(joinpath(data_dir, "small-features/002_scarf.png")))

new = iround(data/255.0)
new = iround(data/255.0 * 0.5); ImageView.display(Image(255 * new))

# ImageView.display(Image(new))

# 4.0 worked kind of
# 3.0 worked quite well
ising = SML.IsingFunction(new, 5.0, 1.0, 1.0, 1.0)
shape = size(new)
A = SML.min_norm_point(ising, [1:prod(shape)], 1e-12)

result = convert(Array{Int, 2}, reshape(A, (shape[2], shape[1])))
img = Image(255*result)
ImageView.display(img)

NPZ.npzwrite("/home/pcmoritz/data/003_comb.npz", result)

# SML.evaluate(ising, SML.ind_to_set(A))
# SML.evaluate(ising, SML.ind_to_set(1.0 * reshape(new, size(A))))

# for file in files  
# end

# python:

import Image
import numpy

d = numpy.load("denoised.npz")
im = Image.fromarray(numpy.uint8((1.-d)*255))
im.convert('RGB').save("denoised.png")

d = numpy.load("noisy.npz")
im = Image.fromarray(numpy.uint8((1.-d)*255))
im.convert('RGB').save("noisy.png")

d = numpy.load("/home/pcmoritz/data/019_redtree1.npz")
im = Image.fromarray(numpy.uint8(d*255))

im.convert('RGB').save("/home/pcmoritz/data/019_redtree1.png")
