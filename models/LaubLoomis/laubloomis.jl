using ReachabilityAnalysis, Plots

# canonical direction along x₄
const e4 = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
const T = 20.0

@taylorize function laubloomis!(dx, x, params, t)
    dx[1] = 1.4*x[3] - 0.9*x[1]
    dx[2] = 2.5*x[5] - 1.5*x[2]
    dx[3] = 0.6*x[7] - 0.8*(x[2]*x[3])
    dx[4] = 2 - 1.3*(x[3]*x[4])
    dx[5] = 0.7*x[1] - (x[4]*x[5])
    dx[6] = 0.3*x[1] - 3.1*x[6]
    dx[7] = 1.8*x[6] - 1.6*(x[2]*x[7])
    return dx
end

function laubloomis(; W=0.01)

    # initial states
    X0c = [1.2, 1.05, 1.5, 2.4, 1.0, 0.1, 0.45]
    X0 = Hyperrectangle(X0c, fill(W, 7))

    # initil-value problem
    prob = @ivp(x' = laubloomis!(x), dim: 7, x(0) ∈ X0)

    return prob
end
