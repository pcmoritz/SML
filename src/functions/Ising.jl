# (c) Philipp Moritz, 2014
# The ising model for image denoising. Let I be a rectangular reference image, then
# F(X) = coeffPix * \sum_{ij} |X_{ij} - I_{ij}| +
#        coeffH * \sum_{ij} \sum_{kl horiz neighbor of ij} |X_{ij} - X_{kl}| +
#        coeffV * \sum_{ij} \sum_{kl vert neighbor of ij} |X_{ij} - X_{kl}| +
#        coeffD * \sum_{ij} \sum_{kl diag neighbor of ij} |X_{ij} - X_{kl}|

type IsingFunction <: Expr
    img::Array{Int, 2}
    mask::Array{Float64, 2}
    height::Int
    width::Int
    value::Float64
    initial_value::Float64
    coeffPix::Float64
    coeffH::Float64
    coeffV::Float64
    coeffD::Float64
end

function IsingFunction(img::Array{Int,2}, coeffPix, coeffH, coeffV, coeffD)
    value = 0.0
    height = Base.size(img, 1)
    width = Base.size(img, 2)
    for j in 1:width
        for i in 1:height
            value += Base.abs(img[i,j])
        end
    end
    return IsingFunction(img, zeros(height,width), height, width, coeffPix*value,
                         coeffPix*value, coeffPix, coeffH, coeffV, coeffD)
end

function reset(func::IsingFunction)
    fill!(func.mask, 0)
    func.value = func.initial_value
end

emptyval(func::IsingFunction) = func.initial_value
currval(func::IsingFunction) = func.value

# Assume the bit at element is not set before this function is called
function incremental(func::IsingFunction, element::Int)
    lastvalue = func.value
    (x, y) = divrem(element-1, func.width)
    x = x + 1
    y = y + 1
    func.mask[x,y] = 1
    if x > 1
        func.value += (1 - 2*func.mask[x-1,y]) * func.coeffH
        if y > 1
            func.value += (1 - 2*func.mask[x-1,y-1]) * func.coeffD
        end
    end
    if  y > 1
        func.value += (1 - 2*func.mask[x, y-1]) * func.coeffV
        if x < func.height
            func.value += (1 - 2*func.mask[x+1,y-1]) * func.coeffD
        end
    end
    if x < func.height
        func.value += (1 - 2*func.mask[x+1,y]) * func.coeffH
        if y < func.width
            func.value += (1 - 2*func.mask[x+1,y+1]) * func.coeffD
        end
    end
    if y < func.width
        func.value += (1 - 2*func.mask[x, y+1]) * func.coeffV
        if x > 1
            func.value += (1 - 2*func.mask[x-1,y+1]) * func.coeffD
        end
    end
    if func.img[x,y] == 1
        func.value -= func.coeffPix
    else
        func.value += func.coeffPix
    end
    return func.value - lastvalue
end

# Andreas Krause's implementation
function evaluate(func::IsingFunction, set::Array{Int})
    image = zeros(Int, func.height, func.width)
    for element in set
        (x, y) = divrem(element-1, func.width)
        x = x + 1
        y = y + 1
        image[x, y] = 1
    end
    # println(image)
    evaluate_image(func, image)
end

function evaluate_image(func::IsingFunction, set::Array{Int, 2})
    delta = Base.abs(set - func.img)
    Epix = sum(delta) # energy through pixelwise disagreement (measured in units of coeffPix)
    # println("Epix")
    # println(Epix)
    Ehor = sum(Base.abs(set[:, 2:func.width] - set[:,1:func.width-1]))
    # println("Ehor")
    # println(Ehor)
    Evert = sum(Base.abs(set[2:func.height,:] - set[1:func.height-1,:]))
    # println("Evert")
    # println(Evert)
    Ediag = sum(Base.abs(set[1:func.height-1, 1:func.width-1] - set[2:func.height, 2:func.width])
        + Base.abs(set[2:func.height, 1:func.width-1] - set[1:func.height-1,2:func.width]))
    # println("Ediag")
    # println(Ediag)
    return func.coeffPix * Epix + func.coeffH * Ehor + func.coeffV * Evert + func.coeffD * Ediag
end
