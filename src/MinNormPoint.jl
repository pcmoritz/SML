# (c) Philipp Moritz, 2014

# The implementation is based on the C implementation by S. Fujishige
# see the following paper:
# S. Fujishige and S. Isotani: A submodular function minimization
# algorithm based on the minimum-norm base. {\it Pacific Journal
# of Optimization} {\bf 7} (2011) 3--17.

const MAX_ITER = 10000

function without_column(A, j)
    return [A[:,1:j-1] A[:,j+1:end]]
end

function without_column!(A, j)
    for k = j:Base.size(A, 2)-1
        for i = 1:Base.size(A,1)
            A[i,k] = A[i,k+1]
        end
    end
    for i = 1:Base.size(A,1)
        A[i,end] = 0.0
    end
end

function extend_column!(A, a, k)
    if Base.size(A,2) >= k
        for i = 1:Base.size(A,1)
            A[i,k] = a[i]
        end
        return A
    end
    return [A a]
end

function tridiagonalize(R, k)
    for j = 1 : k
        a = R[j,j]
        b = R[j+1,j]
        h = sqrt(a*a + b*b)
        c = a/h
        s = b/h
        for jj = j:k
            a = R[j,jj]
            b = R[j+1,jj]
            x = c*a + s*b
            y = c*b - s*a
            R[j,jj] = x
            R[j+1,jj] = y
        end
    end
end

function minimizer(func, init_perm, EPS)
    n = length(init_perm)
    s = zeros(n)
    data = RadixData(n)

    S = zeros(Float64, n, 1)
    R = zeros(Float64, n, 1)

    min_F = polytope_corner!(func, init_perm, s)
    S[:,1] = s
    ss = dot(s, s)
    R[1,1] = sqrt(1.0 + ss)
    k = 1
    kx = 0
    xxp = 1e10
    w = ones(Float64, 1)
    x = zeros(Float64, length(init_perm))
    v = zeros(Float64, 1)
    for i = 1:MAX_ITER
        fill!(x, 0.0)
        for j = 1:k
            for i = 1:n
                x[i] += S[i,j] * w[j]
            end
        end

        min_F = greedy(func, x, data, s)

        ss = dot(s, s)
        xx = dot(x, x)
        sx = dot(s, x)
        
        xxp = xx

        den = max(1.0, ss)

        if (xx - sx)/den < EPS
            return x
        end

        for i = 1:k 
            flag = false
            for j = 1:n
                psik = S[j, i] - S[j, k]
                if psik * psik > EPS
                    flag = true
                    break
                end
            end
            if !flag
                break
            end
        end

        fill!(v, 0.0)        
        for j = 1:k
            for i = 1:n
                v[j] += s[i] * S[i,j]
            end
        end
        for j = 1:k
            v[j] += 1.0
        end
        
        # v = S[:,1:k].' * s +  ones(Float64, k)
        r = LAPACK.trtrs!('U', 'T', 'N', R[1:k,1:k], v)
        # r = R[1:k,1:k].' \ v
        rho = sqrt(1.0 + ss - dot(r, r))
        newcol = [r, rho, zeros(Float64, length(init_perm) - k - 1)]

        # R = [R newcol]
        # S = [S s]
        k += 1
        w = [w, 0.0]
        v = [v, 0.0]
        R = extend_column!(R, newcol, k)
        S = extend_column!(S, s, k)

        # Comment in for debuggin purposes (doesn't work for large matrices)
        # O = ones(Float64, k)
        # assert(sum(abs(O * O.' + S.' * S - R.' * R)) <= EPS)

        while(true)
            rhs = ones(Float64, k)
            # ubar = Triangular(R[1:k,1:k].') \ rhs
            ubar = LAPACK.trtrs!('U', 'T', 'N', R[1:k,1:k], rhs)
            # u = Triangular(R[1:k,1:k]) \ ubar
            u = LAPACK.trtrs!('U', 'N', 'N', R[1:k,1:k], ubar)
            usum = sum(u)
            v = u / usum

            stop_minor = true
            for j = 1:k
                if v[j] < EPS
                    stop_minor = false
                    break
                end
            end

            if stop_minor
                w = v
                break
            end

            direction = w - v
            bounds = w ./ direction

            beta = 1.0
            for j = 1:k
                if direction[j] > EPS
                    beta = min(beta, bounds[j])
                end
            end

            w = w - beta * direction

            j = k
            while j >= 1
                if w[j] < EPS
                    w = [w[1:j-1], w[j+1:end]]
                    v = [v[1:j-1], v[j+1:end]]
                    # history = [history[1:j-1], history[j+1:end]]
                    # S = without_column(S, j)
                    # R = without_column(R, j)
                    without_column!(S, j)
                    without_column!(R, j)
                    tridiagonalize(R, k)
                    k -= 1
                end
                j -= 1
            end
        end
    end
    return x
end

function min_norm_point(func, init_perm, EPS=1e-10)
    p = length(init_perm)
    x = minimizer(func, init_perm, EPS)
    thresh = EPS/(2.0 * p)
    result = zeros(p)
    for i = 1:p
        result[i] = x[i] > -thresh ? 0.0 : 1.0
    end
    return result
end

function minimize(func)
    p = size(func)
    result = min_norm_point(func, [1:p])
    return ind_to_set(result)
end
