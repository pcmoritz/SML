# (c) Philipp Moritz, 2014

# minimize lambda * f(w) + .5 * |w - z|^2
function prox_operator(z, lambda, func; gap = 1e-8)
    f = Prod(Const(lambda), func) + Modular(-z)
    return -minimizer(f, [1:length(z)], gap)
end

# minimize lambda * f(|w|) + .5 * |w - z|^2
# using the proximal operator for f(w) and then project the result
function prox_operator_abs(z, lambda, func; gap=1e-8)
    signs = sign(z)
    z = abs(z)
    w = prox_operator(z, lambda, func; gap=gap)
    return max(w, 0) .* signs
end
