using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
model = "LALO20"
cases = ["W=0.01", "W=0.05", "W=0.1"]
SUITE[model] = BenchmarkGroup()

include("laubloomis.jl")
validation = []
final_width = []

# ----------------------------------------
#  Case 1: smaller initial states
# ----------------------------------------
W = 0.01
prob = laubloomis(W=W)
alg = TMJets(abs_tol=1e-11, orderT=7, orderQ=1, adaptive=true)

# warm-up run
sol_1 = solve(prob, T=T, alg=alg)
sol_1z = overapproximate(sol_1, Zonotope)

# verify that specification holds
property = ρ(e4, sol_1z) < 4.5
push!(validation, Int(property))

# width of final box
width = ρ(e4, sol_1z[end]) + ρ(-e4, sol_1z[end])
push!(final_width, trunc(width, digits=4))
println("width of final box, case $(cases[1]) : $width")

# benchmark
SUITE[model][cases[1]] = @benchmarkable solve($prob, T=$T, alg=$alg)

# ----------------------------------------
# Case 2: intermediate initial states
# ----------------------------------------
W = 0.05
prob = laubloomis(W=W)
alg = TMJets(abs_tol=1e-12, orderT=7, orderQ=1, adaptive=false)

# warm-up run
sol_2 = solve(prob, T=T, alg=alg)
sol_2z = overapproximate(sol_2, Zonotope)

# verify that specification holds
property = ρ(e4, sol_2z) < 4.5
push!(validation, property ? 1 : 0)

# width of final box
width = ρ(e4, sol_2z[end]) + ρ(-e4, sol_2z[end])
push!(final_width, trunc(width, digits=3))
println("width of final box, case $(cases[2]): $width")

# benchmark
SUITE[model][cases[2]] = @benchmarkable solve($prob, T=$T, alg=$alg)

# ----------------------------------------
# Case 3: larger initial states
# ----------------------------------------
W = 0.1
prob = laubloomis(W=W)
alg = TMJets(abs_tol=1e-12, orderT=7, orderQ=1, adaptive=false)

# warm-up run
sol_3 = solve(prob, T=T, alg=alg)
sol_3z = overapproximate(sol_3, Zonotope)

# verify that specification holds
property = ρ(e4, sol_3z) < 5.0
push!(validation, Int(property))

# width of final box
width = ρ(e4, sol_3z[end]) + ρ(-e4, sol_3z[end])
push!(final_width, trunc(width, digits=3))
println("width of final box, case $(cases[3]) : $width")

# benchmark
SUITE[model][cases[3]] = @benchmarkable solve($prob, T=$T, alg=$alg)

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
    print(io, "JuliaReach, $model, $c, $(validation[i]), $(runtimes[c]), $(final_width[i])\n")
end


# ==============================================================================
# Plot
# ==============================================================================

fig = Plots.plot()

Plots.plot!(fig, sol_3z, vars=(0, 4), linecolor="green", color=:green, alpha=0.8)
Plots.plot!(fig, sol_2z, vars=(0, 4), linecolor="blue",  color=:blue, alpha=0.8)
Plots.plot!(fig, sol_1z, vars=(0, 4), linecolor="yellow", color=:yellow, alpha=0.8,
    tickfont=font(30, "Times"), guidefontsize=45,
    xlab=L"t", #\raisebox{-0.5mm}{\textcolor{white}{.}}",
    ylab=L"x_4", # \raisebox{2mm}{\textcolor{white}{.}}",
    xtick=[0., 5., 10., 15., 20.], ytick=[1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5.],
    xlims=(0., 20.), ylims=(1.5, 5.02),
    bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm,
    size=(1000, 1000))

Plots.plot!(fig, x->x, x->4.5, 0., 20., line=2, color="red", linestyle=:dash, legend=nothing)
Plots.plot!(fig, x->x, x->5., 0., 20., line=2, color="red", linestyle=:dash, legend=nothing)

savefig("ARCH-COMP20-JuliaReach-LaubLoomis.png")
