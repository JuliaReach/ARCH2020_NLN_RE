using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
model = "CVDP20"
cases = ["μ=1.0", "μ=2.0"]
SUITE[model] = BenchmarkGroup()

include("vanderpol.jl")
validation = []

# ----------------------------------------
# Case 1: μ = 1
# ----------------------------------------

prob = vanderpolN2(μ=1.0)
alg = TMJets(abs_tol=1e-10, orderT=7, orderQ=1, adaptive=true)

# warm-up run
sol_cvdp1 = solve(prob, T=7.0, alg=alg);
solz_cvdp1 = overapproximate(sol_cvdp1, Zonotope);

# verify that specification holds
property = (ρ(e1y, solz_cvdp1) < 2.75) && (ρ(e2y, solz_cvdp1) < 2.75)
push!(validation, Int(true))

# benchmark
SUITE[model][cases[1]] = @benchmarkable solve($prob, T=7.0, alg=$alg)


# ----------------------------------------
# Case 2: μ = 2
# ----------------------------------------

prob = vanderpolN2(μ=2.0)
alg = TMJets(abs_tol=1e-10, orderT=7, orderQ=1, adaptive=true)

# warm-up run
sol_cvdp2 = solve(prob, T=8.0, alg=alg);
solz_cvdp2 = overapproximate(sol_cvdp2, Zonotope);

# verify that specification holds
property = (ρ(e1y, solz_cvdp2) < 4.05) && (ρ(e2y, solz_cvdp2) < 4.05)
push!(validation, Int(true))

# benchmark
SUITE[model][cases[2]] = @benchmarkable solve($prob, T=8.0, alg=$alg)


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
    print(io, "JuliaReach, $model, $c, $(validation[i]), $(runtimes[c])\n")
end


# ==============================================================================
# Plot
# ==============================================================================

fig = Plots.plot()

fig = plot()
plot!(fig, solz_cvdp1,  vars=(1, 2), lw=0.0, alpha=1.0, color=:red)

plot!(fig, solz_cvdp2,  vars=(1, 2), lw=0.0, alpha=1.0, color=:blue,
    tickfont=font(30, "Times"), guidefontsize=45,
    xlab=L"x_{1}",
    ylab=L"y_1",
    xtick=[-2.0, 0.0, 2.0], ytick=[-4.0, -2.0, 0.0, 2.0, 4.0],
    xlims=(-2.5, 2.5), ylims=(-4.05, 4.05),
    bottom_margin=8mm, left_margin=2mm, right_margin=8mm, top_margin=3mm,
    size=(1000, 1000))

savefig(fig, "ARCH-COMP20-JuliaReach-VanDerPol.png")
