using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
using BenchmarkTools: minimum, median

using ReachabilityAnalysis: BoxEnclosure, BoxClustering, TemplateHullIntersection, LazyClustering

include("spacecraft.jl")

boxdirs = BoxDirections(5)

function solve_spacecraft(prob; k=25, s=missing)

    # transition from mode 1 to mode 2
    sol12 = solve(prob,
                tspan=(0.0, 200.0),
                alg=TMJets(abs_tol=1e-5, max_steps=10_000, orderT=5, orderQ=1, disjointness=BoxEnclosure()),
                max_jumps=1,
                intersect_source_invariant=false,
                intersection_method=TemplateHullIntersection(boxdirs),
                clustering_method=LazyClustering(1),
                disjointness_method=BoxEnclosure())

    sol12jump = overapproximate(sol12[2](120 .. 150), Zonotope)
    t0 = tstart(sol12jump[1])
    sol12jump_c = cluster(sol12jump, 1:length(sol12jump), BoxClustering(k, s))

    # transition from mode 2 to mode 3
    H = system(prob)
    sol3 = solve(IVP(mode(H, 3), [set(X) for X in sol12jump_c]),
                 tspan=(t0, 200.0),
                 alg=TMJets(abs_tol=1e-10, orderT=7, orderQ=1, disjointness=BoxEnclosure()))
    d = Dict{Symbol, Any}(:loc_id => 3)

    return HybridFlowpipe(vcat([fp for fp in sol12.F],
                               [Flowpipe(fp.Xk, d) for fp in sol3.F]))
end

prob = spacecraft()
sol = solve_spacecraft(prob; k=25, s=missing)
solz = overapproximate(sol, Zonotope)

@show line_of_sight(solz)
@show velocity_constraint(solz)
@show target_avoidance(solz)

idx_approaching = findall(x -> x == 1, location.(solz))
idx_attempt = findall(x -> x == 2, location.(solz))
idx_aborting = findall(x -> x == 3, location.(solz))

fig = plot(legend=:bottomright)

for idx in idx_approaching
    plot!(fig, solz[idx], vars=(1, 2), lw=0.0, color=:blue, alpha=1.)
end
for idx in idx_attempt
    plot!(fig, solz[idx], vars=(1, 2), lw=0.0, color=:red, alpha=1.)
end
for idx in idx_aborting
    plot!(fig, solz[idx], vars=(1, 2), lw=0.0, color=:green, alpha=1.)
end
