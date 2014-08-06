# (c) Philipp Moritz, 2014
# This implements the log determinant function; i.e. given a matrix M, 
# the function computes F(X) = \log\det M_X where X is a set of 
# indices that define a submatrix of M.

type LogDet <: Expr
    sigma::Array{Float64, 2}
    indices::BitArray{1}
    chol::Array{Float64, 2}
    det::Float64
    k::Int
end

LogDet(sigma::Array{Float64, 2}) = begin
    n = Base.size(sigma, 1)
    @assert n == Base.size(sigma, 2)
    @assert isposdef(sigma)
    LogDet(sigma, BitArray(n), zeros(n, n), 1.0, 0)
end

reset(expr::LogDet) = begin
    fill!(expr.indices, false)
    fill!(expr.chol, 0.0)
    expr.det = 1.0
    expr.k = 0
end

incremental(expr::LogDet, element::Int) = begin
    expr.k = chol_update!(expr.sigma, expr.indices, element, expr.chol, expr.k)
    expr.indices[element] = true
    return 2.0 * Base.log(expr.chol[expr.k, expr.k])
end

evaluate(expr::LogDet, set::Vector{Int}) = Base.log(det(expr.sigma[set, set]))
size(expr::LogDet) = Base.size(expr.sigma, 1)

emptyval(expr::LogDet) = 0.0

# Based on the implementation of Andreas Krause, originally from Ram
# Rajagopal, see sfo_chol_update.m in the SFO Toolbox
function chol_update!(sigma::Array{Float64, 2}, A, index, chol, k)
    xtx = sigma[index, index]
    rhs = sigma[A, index]
    if k == 0
        k += 1
        chol[k,k] = Base.sqrt(xtx)
        return k
    end
    # Beware: rhs gets changed here
    LAPACK.trtrs!('U', 'T', 'N', chol[1:k,1:k], rhs)
    rpp = xtx - dot(rhs, rhs)
    for j = 1:k
        chol[j,k+1] = rhs[j]
    end
    chol[k+1,k+1] = Base.sqrt(rpp)
    k += 1
    return k
end

# Using tridiagnoalization for cholesky downdate
function chol_downdate!(index, chol, k)
    without_column!(chol, index)
    k -= 1
    tridiagonalize(chol, k)
    return k
end


# assume here that sigma is pos def (see the assert in the constructor)
curvature(expr::LogDet) = :submodular
