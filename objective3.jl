using LinearAlgebra
using Statistics
using Images, FileIO

function gg(f,h,v ; g, λ=1, αw=1)
    """
    geman and geman loss function.
    problems @ edges? I just have component of loss contribute
    zero if it references pixels outside image.
    """
    imheight, imwidth = size(f)
    E = 0
    for i in 1:imheight
        for j in 1:imwidth
            if i != imheight && j != imwidth
                # somewhere in the middle
                E += λ^2 * (f[i+1, j] - f[i,j])^2 *(1 - h[i,j]) +
                λ^2 * (f[i, j+1] - f[i,j])^2*(1 - v[i,j]) +
                αw*(h[i,j] + v[i,j]) +
                (f[i,j] - g[i,j])^2
            elseif i != imheight && j == imwidth
                # right hand side
                E += λ^2 * (f[i+1, j] - f[i,j])^2 *(1 - h[i,j]) +
                αw*h[i,j] +
                (f[i,j] - g[i,j])^2
            elseif i == imheight && j != imwidth
                # bottom row
                E += λ^2 * (f[i, j+1] - f[i,j])^2*(1-v[i,j]) +
                αw*v[i,j] +
                (f[i,j] - g[i,j])^2
            else
                # bottom corner
                E += (f[i,j] - g[i,j])^2
            end
            #print("i = ", i, " j = ", j, "\n")
        end
    end
    return E
end

function get_∇E_wrt_f(f, h, v, λ, image)
    imheight, imwidth = size(f)
    ∇E_wrf_f = zeros(imheight, imwidth)
    for i in 1:imheight
        for j in 1:imwidth
            ∇E_wrf_f[i,j] = get_∇E_wrt_f_single(f, h, v, i, j, λ, image)
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
    ∇E_wrf_v = zeros(h , (w-1))
    for i in 1:h
        for j in 1:(w-1)
            ∇E_wrf_v[i,j] = α - λ^2*(f[i, j+1] - f[i,j])^2
        end
    end
    return ∇E_wrf_v
end



function get_∇E_wrt_h(f, α, λ)
    h, w = size(f)
    ∇E_wrf_h = zeros((h-1) , w)
    for i in 1:(h-1)
        for j in 1:w
            ∇E_wrf_h[i,j] = α - λ^2*(f[i+1, j] - f[i,j])^2
        end
    end
    return ∇E_wrf_h
end


function ∇gg(f,h,v; image, λ=1, αw=1)
    """
    Gradient of geman and geman loss function.
    """
    ∇E_f = get_∇E_wrt_f(f, h, v, λ, image)
    ∇E_v = get_∇E_wrt_v(f, αw, λ)
    ∇E_h = get_∇E_wrt_h(f, αw, λ)
    return ∇E_f, ∇E_h, ∇E_v
end


function p(h,v, ρ)
    """
    This barrier is ρ[xlog(x) +(1-x)log(1-x)]
    Want ρ decreasing from maybe 1 to .1
    """
    hCost = ρ .* sum(h .* log.(h) .+ (1 .- h).*log.(1 .- h) .+ log(2))
    vCost = ρ .* sum(v .* log.(v) .+ (1 .- v).*log.(1 .- v) .+ log(2))
    return vCost + hCost
end

function ∇p(h,v, ρ)
    """
    log(.1) = -2.3, log(.01) = -4.6
    """
    ∇v = ρ .* (log.(v) .- log.(1 .- v))
    ∇h = ρ .* (log.(h) .- log.(1 .- h))
    return ∇h, ∇v
end


function q(h,v, ρ)
    """
    This barrier is ρ[-log(x) - log(1-x)]
    Want ρ decreasing from maybe 1 to .1
    """
    hCost = -ρ .* sum(log.(h) .+ log.(1 .- h) .+ 2*log(2))
    vCost = -ρ .* sum(log.(v) .+ log.(1 .- v) .+ 2*log(2))
    return vCost + hCost
end

function ∇q(h,v, ρ)
    """
    log(.1) = -2.3, log(.01) = -4.6
    """
    ∇v = -ρ .* (1 ./ v .- 1 ./ (1 .- v))
    ∇h = -ρ .* (1 ./ h .- 1 ./ (1 .- h))
    return ∇h, ∇v
end

function get_α_max(v, ∇v, h, ∇h)
    """
    return 1/2 max alpha
    """
    α1h = -(1 .- h) ./ ∇h
    α1v = -(1 .- v) ./ ∇v
    α0h = h ./ ∇h
    α0v = v ./ ∇v
    α = Inf
    for α_mat in [α1h, α1v, α0h, α0v]
        αs = α_mat[α_mat .> 0]
        if length(αs) > 0
            α_min = minimum(αs)
        else
            α_min = Inf
        end

        if α_min < α
            α = α_min
        end
    end
    return α * 0.999
end



function minimize(ggloss, ∇ggloss, b, ∇b, f, h, v, λ, ρ, αw, im)
    for i in 1:2000
        # print("\n\nIteration ", i, "\n")
        # print("Number of h outside bounds = ", sum(h .< 0) + sum(h .> 1), "\n")
        # print("Number of v outside bounds = ", sum(v .< 0) + sum(v .> 1), "\n")
        # print("f(x) = ", gg(f,v,h ; g=im, λ=λ, αw=αw) + q(v,h,ρ), "\n")
        ∇f, ∇h, ∇v = ∇ggloss(f,h,v;image=im, λ=λ, αw=αw)
        ∇h_penalty, ∇v_penalty = ∇b(h,v, ρ)
        ∇v = ∇v .+ ∇v_penalty
        ∇h = ∇h .+ ∇h_penalty
        α = min(get_α_max(v, ∇v, h, ∇h), 1)
        # print("α_max = ", α, "\n")
        objective = α′ -> ggloss(f .- α′*∇f, h .- α′*∇h, v .- α′*∇v ; g=im, λ=λ, αw=αw) + b(h .- α′*∇h, v .- α′*∇v, ρ)
        j = 0
        while j < 5
            if objective(α) < objective(0)
                #print("exited early at j = ", j, "\n")
                break
            else
                α *= .8
                j += 1
            end
        end
        # print("exited at j = ", j, "\n")
        # print("Selected α = ", α, "\n")
        # print("norm(∇f) = ", norm(∇f), "\n")
        # print("norm(∇h) = ", norm(∇h), "\n")
        # print("norm(∇v) = ", norm(∇v), "\n")
        f = f .- α*∇f
        h = h .- α*∇h
        v = v .- α*∇v
    end
    return f, h, v

end


function interior_point(ggloss, ∇ggloss, b, ∇b, f, h, v, λs, ρs, αws, im)
    """
    Need to find good schedule of parameters to go through.
    """
    for j in 1:length(ρs)
        λ = λs[j]
        ρ = ρs[j]
        αw = αws[j]
        f, h, v = minimize(ggloss, ∇ggloss, b, ∇b, f, h, v, λ, ρ, αw, im)
        ∇f, ∇h, ∇v = ∇ggloss(f,h,v; image=im, λ=λ, αw=αw)
        print("\n\nρ = ", ρ, " done\n")
        print("norm(∇f) = ", norm(∇f), "\n")
        print("norm(∇h) = ", norm(∇h), "\n")
        print("norm(∇v) = ", norm(∇v), "\n")
        # cc = 0.9
        # im = cc .* im .+ (1 - cc) .* f
    end
    return f, h, v

end

# ##some basic tests. These numbers seem to decrease so hopefully I did the gradient correctly.

img = load("./im1/1638815998.png")
imgg = Gray.(img);
im = transpose(convert(Array{Float64}, imgg));
imheight, imwidth = size(im)

# ###### Starting Point of optimization
# f = rand(imheight, imwidth);
# h = rand(imheight- 1, imwidth-1);
# v = rand(imheight - 1, imwidth - 1);
f = .5 .+ randn(imheight, imwidth)/10;
h = .5 .+ randn(imheight-1, imwidth)/10;
v = .5 .+ randn(imheight, imwidth-1)/10;

#λ = .5;
# ρ = 2; # also progress to small rho (.001)
#αw = .02; # want large alpha in the beginning, and progress to small alpha
# f,h,v = minimize(gg, ∇gg, q, ∇q, f, h, v, λ, ρ, αw, im)
# λ =1, ρ=.1, αw = 1 sort of works.
# λ =1, ρ=.1, αw = .1 sort of works.
# λ =10, ρ=.1, αw = 10 does not work - at all.
# λ =1, ρ=.01, αw = .1 does not work.
# using q: λ= .5,rho = .001, αw = .04 is good!


# ### Interior Point
niters = 7
ρs = [1.5625*.2^k for k in 2:niters]
αws = repeat([.001], niters-1) # .005 seems to work with p, .002 seems to work well with q.
λs = repeat([.3], niters - 1)
f,h,v = interior_point(gg, ∇gg, q, ∇q, f, h, v, λs, ρs, αws, im)

# ### One round of minimize
# λ = .5
# ρ = .0001
# αw = .005
# f,h,v = minimize(gg, ∇gg, q, ∇q, f, h, v, λ, ρ, αw, im);

edge_image = 1 .- max.(v[1:29, 1:29],h[1:29, 1:29])
edge_image_binary = edge_image .> .5

mosaic(colorview(Gray, im),
colorview(Gray, f),
colorview(Gray, v),
colorview(Gray, h),
colorview(Gray, edge_image),
colorview(Gray, edge_image_binary))
