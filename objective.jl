using LinearAlgebra

g = zeros((5,5))
g[3,3] = 0.5
function gg(x, g, λ=1, α=1)  # Geman Geman energy. Equation 2 in chambolle1995
    I, J = size(g)
    E = 0
    for i in 1:I-1
        for j in 1:J-1
            E += λ^2 * abs(f[i+1,j]-f[i,j])^2 * (1 - v) + λ^2 * abs(f[i,j+1]-f[i,j])^2 * (1 - h) + (h + v)
        end
    end
    return E
end

function gg(x ; g, λ=1, α=1)
    """
    geman and geman loss function.
    problems @ edges? I just have component of loss contribute
    zero if it references pixels outside image.
    """
    imheight, imwidth = size(g)
    # Define f,v,h
    i1 = imwidth*imheight
    i2 = imwidth*imheight + imheight*(imwidth-1)
    f = x[1:i1]
    f = reshape(f, (imheight,imwidth))

    v = x[(i1+1):i2]
    v = reshape(v, (imheight, imwidth-1))

    h = x[(i2+1):end]
    h = reshape(h, (imheight-1, imwidth))

    #print(size(f), size(v), size(h), size(g))
    E = 0
    for i in 1:imheight-1
        for j in 1:imwidth-1
            E += λ^2 * abs(f[i+1,j]-f[i,j])^2 * (1 - v[i, j]) + 
            λ^2 * abs(f[i,j+1]-f[i,j])^2 * (1 - h[i, j]) + 
            α*(h[i, j] + v[i, j]) + 
            (f[i,j] - g[i,j])^2
        end
    end
    return E
end

function get_∇E_wrt_f(f, h, v, λ, g)
    imheight, imwidth = size(f)
    ∇E_wrf_f = zeros(imheight * imheight)
    counter = 1
    for i in 1:imheight
        for j in 1:imwidth
            #print("index i = " ,i,", j = ", j, "\n")
            ∇E_wrf_f[counter] = get_∇E_wrt_f_single(f, h, v, i, j, λ, g)
            counter += 1
        end
    end
    return ∇E_wrf_f
end

function get_∇E_wrt_f_single(f, h, v, i, j, λ, g)
    """
    9 total cases, calculate each separately.
    """
    imheight,imwidth = size(f)
    if i == 1 && j == 1
        return λ^2*(f[i+1, j] - f[i,j])*(-2)*(1 - h[i,j]) +
        λ^2*(f[i, j+1] - f[i,j])*(-2)*(1-v[i,j]) +
        2*(f[i,j] - g[i,j])
    elseif i == 1 && 1 < j < imwidth
        return λ^2*(f[i+1, j] - f[i,j])*(-2)*(1 - h[i,j]) +
        λ^2*(f[i, j+1] - f[i,j])*(-2)*(1-v[i,j]) +
        2*(f[i,j] - g[i,j]) +
        λ^2 *(f[i,j] - f[i, j-1])*(2)*(1 - v[i, j-1])
    elseif i == 1 && j == imwidth
        return λ^2*(f[i+1, j] - f[i,j])*(-2)*(1 - h[i,j]) +
        2*(f[i,j] - g[i,j])
    elseif 1 < i < imheight && j == 1
        return λ^2*(f[i+1, j] - f[i,j])*(-2)*(1 - h[i,j]) +
        λ^2*(f[i, j+1] - f[i,j])*(-2)*(1-v[i,j]) +
        2*(f[i,j] - g[i,j])
    elseif 1 < i < imheight && 1 < j < imwidth
        return λ^2*(f[i+1, j] - f[i,j])*(-2)*(1 - h[i,j]) +
        λ^2 *(f[i,j] - f[i-1, j])*(2)*(1-h[i-1, j]) +
        λ^2 *(f[i, j+1] - f[i,j])*(-2)*(1-v[i,j]) +
        λ^2 * (f[i,j] - f[i, j-1])*(2)*(1-v[i, j-1]) +
        2*(f[i,j] - g[i,j])
    elseif 1 < i < imheight && j == imwidth
        return λ^2*(f[i+1, j] - f[i,j])*(-2)*(1 - h[i,j]) +
        λ^2 *(f[i,j] - f[i-1, j])*(2)*(1-h[i-1, j]) +
        λ^2 * (f[i,j] - f[i, j-1])*(2)*(1-v[i, j-1]) +
        2*(f[i,j] - g[i,j])
    elseif  i == imheight && j == 1
        return λ^2 *(f[i,j] - f[i-1, j])*(2)*(1-h[i-1, j]) +
        λ^2 *(f[i, j+1] - f[i,j])*(-2)*(1-v[i,j]) +
        2*(f[i,j] - g[i,j])
    elseif i == imheight && 1 < j < imwidth
        return λ^2 *(f[i,j] - f[i-1, j])*(2)*(1-h[i-1, j]) +
        λ^2 *(f[i, j+1] - f[i,j])*(-2)*(1-v[i,j]) +
        λ^2 * (f[i,j] - f[i, j-1])*(2)*(1-v[i, j-1]) +
        2*(f[i,j] - g[i,j])
    else # i == h && j == h
        return λ^2 *(f[i,j] - f[i-1, j])*(2)*(1-h[i-1, j]) +
        λ^2 * (f[i,j] - f[i, j-1])*(2)*(1-v[i, j-1]) +
        2*(f[i,j] - g[i,j])
    end
end

function get_∇E_wrt_v(f, α, λ)
    h, w = size(f)
    ∇E_wrf_v = zeros(h * (w-1))
    counter = 1
    for i in 1:h
        for j in 1:(w-1)
            ∇E_wrf_v[counter] = get_∇E_wrt_v_single(f, i, j, α, λ)
            counter += 1
        end
    end
    return ∇E_wrf_v
end

function get_∇E_wrt_v_single(f, i, j, α, λ)
    h, w = size(f)
    if j == w
        return α
    else
        return α - λ^2*(f[i, j+1] - f[i,j])^2
    end
end

function get_∇E_wrt_h(f, α, λ)
    h, w = size(f)
    ∇E_wrf_h = zeros((h-1) * w)
    counter = 1
    for i in 1:(h-1)
        for j in 1:w
            ∇E_wrf_h[counter] = get_∇E_wrt_h_single(f, i, j, α, λ)
            counter += 1
        end
    end
    return ∇E_wrf_h
end

function get_∇E_wrt_h_single(f, i, j, α, λ)
    imheight, w = size(f)
    if i == imheight
        return α
    else
        return α - λ^2*(f[i+1, j] - f[i,j])^2
    end
end


function ∇gg(x; g, λ=1, α=1)
    """
    Gradient of geman and geman loss function.
    """

    imheight, imwidth = size(g)
    # Define f,v,h
    i1 = imwidth*imheight
    i2 = imwidth*imheight + imheight*(imwidth-1)
    f = x[1:i1]
    f = reshape(f, (imheight,imwidth))

    v = x[(i1+1):i2]
    v = reshape(v, (imheight, imwidth-1))

    h = x[(i2+1):end]
    h = reshape(h, (imheight-1, imwidth))
    ∇E_f = get_∇E_wrt_f(f, h, v, λ, g)
    ∇E_v = get_∇E_wrt_v(f, α, λ)
    ∇E_h = get_∇E_wrt_h(f, α, λ)
    ∇E = vcat(∇E_f, ∇E_v, ∇E_h)
    return ∇E
end


function p(x; ρ, g)
    """
    Do I want u nonzero? U for image pixels f, Want to constrain to [0,1],
    but don't want to force towards boundaries like we do with h, v.
    """
    i = prod(size(g))
    u = -1 .* sum(x[1:i] .* log.(x[1:i]) .+ (1 .- x[1:i]).*log.(1 .- x[1:i]))
    #u = 0
    v = -ρ .* sum(x[(i+1):end] .* log.(x[(i+1):end]) .+ (1 .- x[(i+1):end]).*log.(1 .- x[(i+1):end]))
    return u + v
end

function ∇p(x; ρ, g)
    i = prod(size(g))
    u = -1 .* (log.(x[1:i]) .- log.(1 .- x[1:i]))
    #u = zeros(i)
    v = -ρ .* (log.(x[(i+1):end]) .- log.(1 .- x[(i+1):end]))
    return vcat(u,v)
end


function x_to_fhv(x, g)
    imheight, imwidth = size(g)
    # Define f,v,h
    i1 = imwidth*imheight
    i2 = imwidth*imheight + imheight*(imwidth-1)
    f = x[1:i1]
    f = reshape(f, (imheight,imwidth))

    v = x[(i1+1):i2]
    v = reshape(v, (imheight, imwidth-1))

    h = x[(i2+1):end]
    h = reshape(h, (imheight-1, imwidth))
    return (f,h,v)
end

