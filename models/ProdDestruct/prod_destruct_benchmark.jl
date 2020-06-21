using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
model = "PRDE20"
cases = ["I", "P", "IP"]
SUITE[model] = BenchmarkGroup()

include("prod_destruct.jl")
validation = []
final_volume = []

# ----------------------------------------
#  Case 1: uncertain initial conditions
# ----------------------------------------

prob = production_destruction(case=1)
alg = TMJets(abs_tol=1e-11, orderT=7, orderQ=1, adaptive=true)

# warm-up run
sol_pd1 = solve(prob, T=100.0, alg=alg)
sol_pd1z = overapproximate(sol_pd1, Zonotope)

# verify that specification holds
property, vol = prod_dest_property(sol_pd1z)
push!(validation, Int(property))
push!(final_volume, trunc(vol, sigdigits=2))

println("volume final box, case $(cases[1]) : $(final_volume[1])")

# benchmark
SUITE[model][cases[1]] = @benchmarkable solve($prob, T=100.0, alg=$alg)


# ----------------------------------------
#  Case 2: uncertain parameter
# ----------------------------------------

prob = production_destruction(case=2)
alg = TMJets(abs_tol=1e-12, orderT=7, orderQ=1, adaptive=true)

# warm-up run
sol_pd2 = solve(prob, T=100.0, alg=alg)
sol_pd2z = overapproximate(sol_pd2, Zonotope)

# verify that specification holds
property, vol = prod_dest_property(sol_pd2z)
push!(validation, Int(property))
push!(final_volume, trunc(vol, sigdigits=2))

println("volume final box, case $(cases[2]) : $(final_volume[2])")

# benchmark
SUITE[model][cases[2]] = @benchmarkable solve($prob, T=100.0, alg=$alg)



# ----------------------------------------
# Case 3: uncertain initial states and parameter
# ----------------------------------------
prob = production_destruction(case=3)
alg = TMJets(abs_tol=1e-11, orderT=7, orderQ=1, adaptive=true)

# warm-up run
sol_pd3 = solve(prob, T=100.0, alg=alg)
sol_pd3z = overapproximate(sol_pd3, Zonotope)

# verify that specification holds
property, vol = prod_dest_property(sol_pd3z)
push!(validation, Int(property))
push!(final_volume, trunc(vol, sigdigits=2))

println("volume final box, case $(cases[3]) : $(final_volume[3])")

# benchmark
SUITE[model][cases[3]] = @benchmarkable solve($prob, T=100.0, alg=$alg)


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
    print(io, "JuliaReach, $model, $c, $(validation[i]), $(runtimes[c]), $(final_volume[i])\n")
end

#
# ==============================================================================
# Plot
# ==============================================================================

fig = Plots.plot()

dt = 0 .. 100

Plots.plot!(fig, sol_pd3z(dt),  vars=(0, 3), linecolor="red", color=:red, alpha=3.0, legend=:bottomright,
            lab="I & P")

Plots.plot!(fig, sol_pd2z(dt), vars=(0, 3), linecolor="blue", color=:blue, alpha=0.8,
    tickfont=font(30, "Times"), guidefontsize=45,
    xlab=L"t",
    ylab=L"z",
    xtick=[0., 25., 50., 75., 100.], ytick=[0.0, 2.5, 5.0, 7.5, 10.0],
    xlims=(0., 100.5), ylims=(0.0, 11.0),
    bottom_margin=6mm, left_margin=2mm, right_margin=6mm, top_margin=3mm,
    size=(1000, 1000), lab="P", legendfontsize=20)

Plots.plot!(fig, sol_pd1z(dt), vars=(0, 3), linecolor="yellow", color=:yellow, alpha=0.3, lab="I")

savefig("ARCH-COMP20-JuliaReach-ProductionDestruction.png")
