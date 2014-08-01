# Encoding of the Disciplined Submodular Ruleset
# (c) Philipp Moritz, 2014
#
# for the Disciplined Convex Ruleset, see Chapter in Global
# Optimization: From Theory to Implementation, L. Liberti and
# N. Maculan (eds.), in the book series Nonconvex Optimization and its
# Applications, Springer, 2006, pages
# 155-210. (http://stanford.edu/~boyd/papers/pdf/disc_cvx_prog.pdf)

using Match

# THE CURVATURE RULES

# Add two curvatures
function plus(first, second)
    @match (first, second) begin
        (:modular, second) => second
        (first, :modular) => first
        (:submodular, :submodular) => :submodular
        (:supmodular, :supmodular) => :supmodular
        (_, _) => :unknown
    end
end

# Negate a curvature
function negative(curvature)
    @match curvature begin
        :submodular => :supmodular
        :supmodular => :submodular
        :modular => :modular
        :unknown => :unknown
    end
end

# Multiply curvature and scalar
function multiply(scalar, curvature)
    if scalar < 0.0
        return negative(curvature)
    else
        return curvature
    end
end

function non_decreasing(mon)
    return mon == :constant || mon == :increasing
end

function non_increasing(mon)
    return mon == :constant || mon == :decreasing
end

# curvatures could be :submodular, :supmodular or :modular
# monotinicities could be :increasing, :decreasing or :constant
function composition(outer_mon::Vector{Symbol}, inner_curv::Vector{Symbol}, inner_mon::Vector{Symbol})
    @assert length(inner_curv) == length(inner_mon)
    @assert length(outer_mon) == length(inner_curv)

    result_submodular = true
    
    for i = 1:length(outer_mon)
        if inner_curv[i] == :modular && non_decreasing(inner_mon[i])
            continue
        end
        if non_decreasing(outer_mon[i])
            result_submodular &= non_decreasing(inner_mon[i]) && inner_curv[i] == :submodular
        end
        if non_increasing(outer_mon[i])
            result_submodular &= non_increasing(inner_mon[i]) && inner_curv[i] == :supmodular
        end
    end

    return result_submodular
end
