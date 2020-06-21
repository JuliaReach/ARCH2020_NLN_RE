# =================================================================
# Quadrotor model
# See https://easychair.org/publications/paper/gjfh
# =================================================================

using ReachabilityAnalysis, Plots

const RA = ReachabilityAnalysis
const T_lv = 3.64
const B = Ball2([1.0, 1.0], 0.15) # "exact"


@taylorize function lotka_volterra!(du, u, p, t)
    u1u2 = u[1] * u[2]
    du[1] = 3.0 * (u[1] - u1u2)
    du[2] = u1u2 - u[2]
    return du
end


function lotka_volterra_hybrid(; nsplit=1,
                                 ε = 0.008,
                                 ε_ext=1e-4, # threshold for the outer approximation
                                 n_int=50)   # number of directions for the inner approximation

    # generate external / internal polytopic approximations of the guard
    B_ext = overapproximate(B, ε_ext) # outer approximation
    B_int = underapproximate(B, PolarDirections(n_int)) # inner approximation
    B_int = tohrep(convert(VPolygon, B_int)) # cast to Hrep
    B_intᶜ = complement(B_int)

    # define modes
    aut = LightAutomaton(3)
    outside = @system(x' = lotka_volterra!(x), dim: 2, x ∈ B_intᶜ)
    inside = @system(x' = lotka_volterra!(x), dim: 2, x ∈ B_ext)
    outside_unconstrained = @system(x' = lotka_volterra!(x), dim: 2, x ∈ Universe(2))

    # define the transition graph
    add_transition!(aut, 1, 2, 1)
    add_transition!(aut, 2, 3, 2)
    T_out_in = @map(x -> x, dim:2, x ∈ B_ext)
    T_in_out = @map(x -> x, dim:2, x ∈ B_intᶜ)

    # initial-value problem
    H = HybridSystem(automaton=aut, modes=[outside, inside, outside_unconstrained],
                                           resetmaps=[T_out_in, T_in_out])

    # initial states with splitting
    X0 = Hyperrectangle(low=[1.3-ε, 1.], high=[1.3+ε, 1.])
    X0s = split(X0, [nsplit, 1])
    X0st = [(X0s_i, 1) for X0s_i in X0s]

    return InitialValueProblem(H, X0st)
end


function lotka_volterra()
    # hybrid problem
    prob = lotka_volterra_hybrid(nsplit=5, ε_ext=1e-6, n_int = 50, ε = 0.008);

    return prob
end


@inline function lv_property(solz; ε_ext=1e-4)

    # Sets intersecting the nonlinear guard
    B = Ball2([1.0, 1.0], 0.15) # "exact"
    B_ext = overapproximate(B, ε_ext) # outer approximation
    intersecting_reachsets = []
    for (i, Fi) in enumerate(solz)
        for (j, Xj) in enumerate(Fi)
            !isdisjoint(Xj, B_ext) && push!(intersecting_reachsets, (i, j))
        end
    end

    # Compute time spent inside non-linear guard
    times = [tspan(solz[ind[1]][ind[2]]) for ind in intersecting_reachsets];
    tmin = minimum(tstart, times)
    tmax = maximum(tend, times)
    @show(tmin, tmax, tmax-tmin)

    indxs = Int[]
    for (i, e) in enumerate(tspan.(solz_lv))
        T_lv ∈ tspan(e) && push!(indxs, i)
    end
    chlast = ConvexHullArray([set(solz[i](T_lv)) for i in indxs]);
    chlasth = overapproximate(chlast, Hyperrectangle)
    a = low(chlasth)
    b = high(chlasth)
    return (b[1] - a[1]) * (b[2] - a[2]), tmax-tmin
end
