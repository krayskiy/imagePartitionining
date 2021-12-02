function interior_point_method(f, p, x; ρ=1, γ=2, ε=0.001)
    delta = Inf
    while delta > ε
        x′ = minimize(x -> f(x) + p(x)/ρ, x)
        delta = norm(x′ - x)
        x = x′
        ρ *= γ
    end
    return x
end