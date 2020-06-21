using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
model = "LOVO20"
cases = [""]
SUITE[model] = BenchmarkGroup()

include("lotka_volterra.jl")
validation = []
final_area = []
intersect_time = []

# ----------------------------------------
#  Case 1:
# ----------------------------------------
prob = lotka_volterra()
alg = TMJets(abs_tol=1e-14, orderT=7, orderQ=1, adaptive=true,
            disjointness=RA.ZonotopeEnclosure())

# warm-up run
sol_lv = solve(prob, T=T_lv,
                  alg=alg,
                  max_jumps=2,
                  intersect_source_invariant=false,
                  intersection_method=RA.BoxIntersection(),
                  clustering_method=RA.BoxClustering(),
                  disjointness_method=RA.BoxEnclosure());
solz_lv = overapproximate(sol_lv, Zonotope);

# obtain area
area, time_in_guard = lv_property(solz_lv)
push!(validation, Int(true))
push!(final_area, trunc(area, sigdigits=3))
push!(intersect_time, trunc(time_in_guard, sigdigits=3))
println("Final area, case $(cases[1]) : $(area)")
println("Time spent in guard, case $(cases[1]) : $(time_in_guard)")

# benchmark
SUITE[model][cases[1]] = @benchmarkable solve($prob,
                  T = $T_lv,
                  alg = $alg,
                  max_jumps = 2,
                  intersect_source_invariant = false,
                  intersection_method = RA.BoxIntersection(),
                  clustering_method = RA.BoxClustering(),
                  disjointness_method = RA.BoxEnclosure())


# ==============================================================================
# Execute benchmarks and save benchmark results
# ==============================================================================

# tune parameters
tune!(SUITE)

# run the benchmarks
results = run(SUITE, verbose=true)

# return the sample with the smallest time value in each test
println("minimum time for each benchmark:\n", minimum(results))

# return the median for each test
println("median time for each benchmark:\n", median(results))

# export runtimes
runtimes = Dict()
for (i, c) in enumerate(cases)
    t = median(results[model][c]).time * 1e-9
    runtimes[c] = t
end

for (i, c) in enumerate(cases)
    print(io, "JuliaReach, $model, $c, $(validation[i]), $(runtimes[c])," *
        " $(final_area[i]), $(intersect_time[i])\n")
end


# ==============================================================================
# Plot
# ==============================================================================

fig = Plots.plot()

fig = plot()
plot!(fig, solz_lv[1],  vars=(1, 2), lw=0.0, alpha=1.0, ratio=1.0, color=:blue)
plot!(fig, solz_lv[2],  vars=(1, 2), lw=0.0, alpha=1.0, ratio=1.0, color=:blue)
plot!(fig, solz_lv[3],  vars=(1, 2), lw=0.0, alpha=1.0, ratio=1.0, color=:blue)
plot!(fig, solz_lv[6],  vars=(1, 2), lw=0.0, alpha=1.0, ratio=1.0, color=:blue)

plot!(fig, solz_lv[4],  vars=(1, 2), lw=0.0, alpha=1.0, ratio=1.0, color=:lightgreen)
plot!(fig, solz_lv[5],  vars=(1, 2), lw=0.0, alpha=1.0, ratio=1.0, color=:lightgreen)
plot!(fig, solz_lv[7],  vars=(1, 2), lw=0.0, alpha=1.0, ratio=1.0, color=:lightgreen)

plot!(fig, B, 1e-4, color=:white, lw=2.0, linecolor=:red, tickfont=font(30, "Times"), guidefontsize=45,
    xlab=L"x",
    ylab=L"y",
    xtick=[0.8, 1.0, 1.2, 1.4], ytick=[0.6, 0.8, 1.0, 1.2, 1.4],
    xlims=(0.6, 1.4), ylims=(0.6, 1.4),
    bottom_margin=6mm, left_margin=2mm, right_margin=8mm, top_margin=3mm,
    size=(1000, 1000))

savefig("ARCH-COMP20-JuliaReach-Lotka-Volterra.png")
