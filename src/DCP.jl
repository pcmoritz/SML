# Encoding of the Disciplined Convex Programming Ruleset
# (c) Philipp Moritz, 2014
#
# see Chapter in Global Optimization: From Theory to Implementation,
# L. Liberti and N. Maculan (eds.), in the book series Nonconvex
# Optimization and its Applications, Springer, 2006, pages
# 155-210. (http://stanford.edu/~boyd/papers/pdf/disc_cvx_prog.pdf)

using Match

# THE CURVATURE RULES

# Add two curvatures
function cvx_plus(first, second)
    @match (first, second) begin
        (:affine, second) => second
        (first, :affine) => first
        (:convex, :convex) => :convex
        (:concave, :concave) => :concave
        (_, _) => :unknown
    end
end

# Add two monotonicities
function mon_plus(first, second)
    @match (first, second) begin
        (:constant, second) => second
        (first, :constant) => first
        (:increasing, :increasing) => :increasing
        (:decreasing, :decreasing) => :decreasing
    end
end

# Negate a curvature
function cvx_negative(curvature)
    @match curvature begin
        :convex => :concave
        :concave => :convex
        :affine => :affine
        :unknown => :unknown
    end
end

# Negate a monotonicity
function mon_negative(monotonicity)
    @match monotonicity begin
        :constant => :constant
        :increasing => :decreasing
        :decreasing => :increasing
    end
end

# Multiply curvature and scalar
function cvx_multiply(scalar, curvature)
        if scalar < 0.0
            return cvx_negative(curvature)
        else
            return curvature
        end
end

# Multiply monotonicity and scalar
function mon_multiply(scalar, monotonicity)
    if scalar < 0.0
        return mon_negative(monotonicity)
    else
        return monotonicity
    end
end

