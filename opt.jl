# Optimizer julia file
include("objective.jl")


function get_α(x,d,f)
    """
    Constrained Line Search
    want 0 < x_i < 1 for all x_i
    Therefore want 0 < x_i + α*d < 1
    Find maximum α such that this not exceeded for any x_i.
    Also, my selection of α is based on grid search rather than anything sophisticated.
    """
    αs = zeros(length(x))
    for i in 1:length(x)
        if d[i] > 0
            # x would increase for positive α, therefore max out at 1
            αs[i] = (1 - x[i])/d[i]
        elseif d[i] < 0
            # x would decrease for positive α, therefore check at zero.
            αs[i] = -x[i]/d[i]
        else
            αs[i] = 1
        end
    end
    max_α = minimum(αs)

    n = 30
    best_α = max_α
    best_f = Inf
    L = LinRange(0,max_α, n)
    for j in 1:n
        fx = f(x .+ L[j]*d)
        if fx < best_f
            best_α = L[j]
            best_f = fx
        end
    end
    return best_α, n

end

# The following code is from K&W, page 74
abstract type DescentMethod end

mutable struct ConjugateGradientDescent <: DescentMethod
    g
    d
end

function init!(M::ConjugateGradientDescent, f, ∇f, x)
    M.g = ∇f(x)
    M.d = -M.g
    return M
end

function step!(M::ConjugateGradientDescent, f, ∇f, x)
    d,g = M.d, M.g
    g2 = ∇f(x)
    β = max(0, dot(g2, g2 .- g)/(g⋅g))
    d2 = -g2 .+ β*d
    #x2, fev = backtracking_line_search_cg(f,g2,x,d2)
    α, fev = get_α(x,d2, f)
    x2 = x .+ α*d2
    M.d, M.g = d2, g2
    return x2, fev+1
end


function minimize(f, ∇f, x, eps=1e-5)
    M = ConjugateGradientDescent(f(x),-∇f(x))
    M = init!(M, f, ∇f, x)
    fev = 0
    convergence = 1
    while convergence > eps
        xnext, fev2 = step!(M, f, ∇f, x)
        convergence = norm(xnext - x)
        x = xnext
        fev += fev2
        #print("\nf(x) = ", f(x), "\n")
    end

    return x, fev
end


function interior_point_method(f, ∇f, p, ∇p, x, image; ρ=1, γ=2, ϵ=1e-6)
    delta = Inf
    while delta > ϵ
        objective = x -> f(x, g=image, λ=1, α=1) + p(x, ρ=ρ, g=image)
        grad = x -> ∇f(x, g=image, λ=1, α=1) .+ ∇p(x, ρ=ρ, g=image)
        x′, fev = minimize(objective, grad, x)
        delta = norm(x′ - x)
        x = x′
        ρ *= γ
    end
    return x
end

# x = .5*ones(280);
# image = rand(10,10)/10;
# image[3:5, 4:7] .= .75;
# f = rand(10,10);
# h = rand(9, 10);
# v = rand(10,9);
# ρ=1;
# objective = x -> gg(x, g=image, λ=1, α=1) + p(x, ρ=ρ, g=image)
# grad = x -> ∇gg(x, g=image, λ=1, α=1) .+ ∇p(x, ρ=ρ, g=image)
# #minx, fev = minimize(objective, grad, x);
# xip = interior_point_method(gg, ∇gg, p, ∇p, x, image);
# # print("minimum(final x) = ", minimum(minx), "\n")
# # print("maximum(final x) = ", maximum(minx), "\n")
# # print("\noriginal Function evaluation = ", objective(x), "\n")
# # print("\nFinal Function evaluation = ", objective(minx), "\n")
# f,h,v = x_to_fhv(xip, image)
