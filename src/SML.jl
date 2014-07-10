# (c) Philipp Moritz, 2014

module SML

include("Utils.jl")
include("DSR.jl")
include("Expr.jl")

include("functions/Modular.jl")
include("functions/Iwata.jl")
include("functions/GraphCut.jl")
include("functions/GraphIntensity.jl")
include("functions/LogDet.jl")
include("functions/SetCover.jl")
include("functions/TreeCover.jl")
include("functions/Ising.jl")
include("functions/Composition.jl")
include("functions/Induced.jl")

include("extras/ProxOperator.jl")

include("Sort.jl")
include("Greedy.jl")
include("MinNormPoint.jl")
# include("Queyranne.jl")

end
