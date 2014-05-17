# (c) Philipp Moritz, 2014

using SML
using PyPlot

maxiter = 100

d = 7
p = 2^d - 1
n = 100
noise = 0.1

function wavelets(d, x)
    X = zeros(Float64, length(x), 2^d-1)
    k = 1
    for i = 1:d
        for j = 1:2:2^i-1
            pos = map(x -> x > (j-1)/2^i && (x <= j/2^i), x)
            neg = map(x -> x > j/2^i && x <= (j+1)/2^i, x)
            X[:,k] = (pos - neg) * 2^(i/2)
            k += 1
        end
    end
    return X
end

F = SML.TreeCoverFunction(d)

x = [0:(n-1)] / (n-1)
X = wavelets(d, x)

plot(x, sin(8*pi*x.^2), "-r")

y = sin(8*pi*x.^2) + noise * randn(n)
ytest = sin(8*pi*x.^2) + noise * randn(n)
meany = mean(y)
yc = y .- meany

plot(x, yc, ".")

L = eigmax(X'*X/n)

function regression(X, y, L, lambda, EPS=1e-12)
    (n, p) = Base.size(X)
    w = zeros(p)
    v = w
    t = 1
    for iter = 1:maxiter
        gradient = 1/n * X' * (X * v - y)
        wold = w
        w = SML.prox_operator_abs(v-1/L*gradient,lambda/L,F;gap=1e-1 * EPS)
        delta = norm(w - wold)
        print("iteration ")
        print(iter)
        print(" ")
        println(delta)
        if(delta < EPS)
            break
        end
        told = t
        t = (1 + sqrt(1 + 4.0*t*t ))/2.0
        v = w + (told - 1) / t * (w - wold)
    end
    return (norm(ytest .- X*w .- meany), w) # return test error
end

lambdas = 10.^[1:-.5:-4.5]

min_loss = typemax(Float64)
best_lambda = typemax(Float64)

for lambda in lambdas
    println("lambda ", lambda)
    (loss, w) = regression(X, yc, L, lambda)
    if loss < min_loss
        min_loss = loss
        best_lambda = lambda
    end
end

(loss, w) = regression(X, yc, L, best_lambda)

plot(x, sin(8*pi*x.^2), "-r")
plot(x, yc, ".")
plot(x, X * w .+ mean(y), "-")
