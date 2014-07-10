using Images
using ImageView
using SML

data_dir = ENV["DATA_DIR"]
directory = joinpath(data_dir, "small-features")
files = readdir(directory)

data = convert(Array, imread(joinpath(directory, files[4])))
ImageView.display(Image(data))

new = iround(data/255.0)

ising = SML.IsingFunction(new, 5.0, 1.0, 1.0, 1.0)
shape = size(new)
A = SML.min_norm_point(ising, [1:prod(shape)], 1e-12)

result = convert(Array{Int, 2}, reshape(A, (shape[2], shape[1])))
img = Image(255*result)
ImageView.display(img)

# SML.evaluate(ising, SML.ind_to_set(A))
# SML.evaluate(ising, SML.ind_to_set(1.0 * reshape(new, size(A))))

# for file in files  
# end
