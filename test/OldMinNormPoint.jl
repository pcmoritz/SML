# (c) Philipp Moritz, 2014

# The implementation is based on the C implementation by S. Fujishige
# that is described in the following paper:
# S. Fujishige and S. Isotani: A submodular function minimization
# algorithm based on the minimum-norm base. {\it Pacific Journal
# of Optimization} {\bf 7} (2011) 3--17.

function without_column(A, j)
    return [A[:,1:j-1] A[:,j+1:end]]
end

function tridiagonalize(R)
    for j = 1 : Base.size(R,2)
        a = R[j,j]
        b = R[j+1,j]
        h = sqrt(a*a + b*b)
        c = a/h
        s = b/h
        for jj = j:Base.size(R,2)
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

    S = zeros(Float64, length(init_perm), 1)
    R = zeros(Float64, length(init_perm), 1)

    decrease = []

    (s, min_F) = polytope_corner(func, init_perm)

    S[:,1] = s
    ss = dot(s, s)

    R[1,1] = sqrt(1.0 + ss)
    k = 1
    kx = 0
    xxp = 1e10
    w = ones(Float64, 1)
    x = zeros(Float64, length(init_perm))

    while(true)
        # x = S * w
        A_mul_B!(x, S, w)
        
        (s, min_F) = greedy(func, x)

        ss = dot(s, s)
        xx = dot(x, x)
        sx = dot(s, x)

        xxp = xx

        den = max(1.0, ss)

        if (xx - sx)/den < EPS
            if xx - sx < 0.0
                error("Are you sure your function is submodular?")
            end
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
        
        v = S.' * s +  ones(Float64, k)
        r = LAPACK.trtrs!('U', 'T', 'N', R[1:k,1:k], v)
        # r = R[1:k,1:k].' \ v
        rho = sqrt(1.0 + ss - dot(r, r))
        newcol = [r, rho, zeros(Float64, length(init_perm) - k - 1)]

        R = [R newcol]
        S = [S s]
        w = [w, 0.0]

        k += 1

        # Comment in for debuggin purposes (doesn't work for large matrices)
        ## O = ones(Float64, k)
        ## if sum(abs(O * O.' + S.' * S - R.' * R)) >= EPS
        ##     println(O * O.' + S.' * S - R.' * R)
        ## end
        ## assert(sum(abs(O * O.' + S.' * S - R.' * R)) <= EPS)
        

        while(true)
            rhs = ones(Float64, k)
            # ubar = Triangular(R[1:k,1:k].') \ rhs
            ubar = LAPACK.trtrs!('U', 'T', 'N', R[1:k,1:k], rhs)
            # u = Triangular(R[1:k,1:k]) \ ubar
            u = LAPACK.trtrs!('U', 'N', 'N', R[1:k,1:k], rhs)
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
                decrease = [decrease; norm(x)]
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
                    S = without_column(S, j)
                    R = without_column(R, j)
                    tridiagonalize(R)
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
