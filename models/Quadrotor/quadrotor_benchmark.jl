using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
model = "QUAD20"
cases = ["Δ=0.1", "Δ=0.4", "Δ=0.8"]
SUITE[model] = BenchmarkGroup()

include("quadrotor.jl")
validation = []

# ----------------------------------------
#  Case 1: smaller uncertainty
# ----------------------------------------
Wpos = 0.1
Wvel = 0.1
prob = quadrotor(project_reachset=false, Wpos=Wpos, Wvel=Wvel)
alg = TMJets(abs_tol=1e-7, orderT=5, orderQ=1, adaptive=false)

# warm-up run
sol1 = solve(prob, tspan=Tspan, alg=alg);
solz1 = overapproximate(sol1, Zonotope);

# verify that specification holds
property = quad_property(solz1)
push!(validation, Int(property))

println("Validate property, case $(cases[1]) : $(property)")

# benchmark
SUITE[model][cases[1]] = @benchmarkable solve($prob, T=$Tspan, alg=$alg)


# ----------------------------------------
# Case 2: intermediate uncertainty
# ----------------------------------------
Wpos = 0.4
Wvel = 0.4
prob = quadrotor(project_reachset=false, Wpos=Wpos, Wvel=Wvel)
alg = TMJets(abs_tol=1e-7, orderT=5, orderQ=1, adaptive=false)

# warm-up run
sol2 = solve(prob, tspan=Tspan, alg=alg);
solz2 = overapproximate(sol2, Zonotope);

# verify that specification holds
property = quad_property(solz2)
push!(validation, Int(property))

println("Validate property, case $(cases[2]) : $(property)")

# benchmark
SUITE[model][cases[2]] = @benchmarkable solve($prob, T=$Tspan, alg=$alg)

# ----------------------------------------
# Case 3: large uncertainty
# ----------------------------------------
Wpos = 0.8
Wvel = 0.8
prob = quadrotor(project_reachset=false, Wpos=Wpos, Wvel=Wvel)
alg = TMJets(abs_tol=1e-7, orderT=5, orderQ=1, adaptive=false)

# warm-up run
sol3 = solve(prob, tspan=Tspan, alg=alg);
solz3 = overapproximate(sol3, Zonotope);

# verify that specification holds
property = quad_property(solz3)
push!(validation, Int(property))

println("Validate property, case $(cases[3]) : $(property)")

# benchmark
SUITE[model][cases[3]] = @benchmarkable solve($prob, T=$Tspan, alg=$alg)


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

Plots.plot!(solz3,  vars=(0, 3), linecolor="green", color=:green, alpha=0.8)
Plots.plot!(solz2, vars=(0, 3), linecolor="blue",   color=:blue, alpha=0.8)
Plots.plot!(solz1, vars=(0, 3), linecolor="yellow",  color=:yellow, alpha=0.8,
    tickfont=font(30, "Times"), guidefontsize=45,
    xlab=L"t", # \raisebox{-0.5mm}{\textcolor{white}{.}}",
    ylab=L"x_3", # \raisebox{2mm}{\textcolor{white}{.}}",
    xtick=[0., 1., 2., 3., 4., 5.], ytick=[-1., -0.5, 0., 0.5, 1., 1.5],
    xlims=(0., 5.), ylims=(-1., 1.5),
    bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm,
    size=(1000, 1000))

Plots.plot!(x->x, x->0.98, 0., 5., line=2, color="red", linestyle=:dash, legend=nothing)
Plots.plot!(x->x, x->1.02, 0., 5., line=2, color="red", linestyle=:dash, legend=nothing)
Plots.plot!(x->x, x->0.9,  0., 5., line=2, color="red", linestyle=:dash, legend=nothing)
Plots.plot!(x->x, x->1.4,  0., 5., line=2, color="red", linestyle=:dash, legend=nothing)

savefig("ARCH-COMP20-JuliaReach-Quadrotor.png")
