# =================================================================
# Coupled Van der Pol model
# See https://easychair.org/publications/paper/nrdD
# =================================================================

using ReachabilityAnalysis, Plots

const e1y = [0.0, 1.0, 0.0, 0.0]
const e2y = [0.0, 0.0, 0.0, 1.0]

@taylorize function vdp_N2_mu1!(dx, x, params, t)
    local μ = 1.0
    x₁, y₁, x₂, y₂ = x

    aux0 = x₂ - x₁
    #
    aux11 = 1 - x₁^2
    aux12 = μ * y₁
    aux13 = aux11 * aux12
    aux15 = aux13 + aux0
    #
    aux21 = 1 - x₂^2
    aux22 = μ * y₂
    aux23 = aux21 * aux22
    aux25 = aux23 - aux0

    dx[1] = y₁
    dx[2] = aux15 - x₁ # (1 - x₁^2)*(μ * y₁) + x₂ - x₁ - x₁
    dx[3] = y₂
    dx[4] = aux25 - x₂ # (1 - x₂^2) * (μ * y₂) - x₂ + x₁ - x₂

    return dx
end

@taylorize function vdp_N2_mu2!(dx, x, params, t)
    local μ = 2.0
    x₁, y₁, x₂, y₂ = x

    aux0 = x₂ - x₁
    #
    aux11 = 1 - x₁^2
    aux12 = μ * y₁
    aux13 = aux11 * aux12
    aux15 = aux13 + aux0
    #
    aux21 = 1 - x₂^2
    aux22 = μ * y₂
    aux23 = aux21 * aux22
    aux25 = aux23 - aux0

    dx[1] = y₁
    dx[2] = aux15 - x₁ # (1 - x₁^2)*(μ * y₁) + x₂ - x₁ - x₁
    dx[3] = y₂
    dx[4] = aux25 - x₂ # (1 - x₂^2) * (μ * y₂) - x₂ + x₁ - x₂

    return dx
end


function vanderpolN2(; μ=1.0)

    if μ == 1.0
        # Initial conditions
        X0 = convert(Hyperrectangle, (1.25..1.55) × (2.35..2.45) × (1.55..1.85) × (2.35..2.45))
        X0 = Hyperrectangle(Vector(X0.center), Vector(X0.radius))
        # initial-value problem
        prob = @ivp(x' = vdp_N2_mu1!(x), dim: 4, x(0) ∈ X0)
    elseif μ == 2.0
        # Initial conditions
        X0 = convert(Hyperrectangle, (1.55..1.85) × (2.35..2.45) × (1.55..1.85) × (2.35..2.45))
        X0 = Hyperrectangle(Vector(X0.center), Vector(X0.radius))
        X0 = split(X0, [7, 1, 7, 1])
        # initial-value problem
        prob = @ivp(x' = vdp_N2_mu2!(x), dim: 4, x(0) ∈ X0)
    else
        error("the value of μ = $μ is not implemented")
    end

    return prob
end
