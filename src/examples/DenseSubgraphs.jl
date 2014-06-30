using MAT
using Graphs

# Download the data file

cd(Pkg.dir("SML") * "/data/")

if !isfile("cnr-2000.mat")
    Base.run(`wget http://www.cise.ufl.edu/research/sparse/mat/LAW/cnr-2000.mat`)
end

# Get the sparse data matrix

vars = MAT.matread("cnr-2000.mat")
A = vars["Problem"]["A"]

# Construct the connectivity graph

g = SML.from_weight_matrix(-A; is_directed=false)

intensity = SML.GraphIntensity(g)

SML.min_norm_point(intensity, [1:num_vertices(g)], 1e-12)




