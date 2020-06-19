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

#=
# ==============================================================================
# Create plots
# ==============================================================================

# --------------------------
# Case 1
# --------------------------

plot(sol_1,
     tickfont=font(30, "Times"), guidefontsize=45,
     xlab=L"x\raisebox{-0.5mm}{\textcolor{white}{.}}",
     ylab=L"y\raisebox{2mm}{\textcolor{white}{.}}",
     xtick=[-3., -2., -1., 0., 1., 2., 3.], ytick=[-3., -2., -1., 0., 1., 2., 3.],
     xlims=(-3., 3.), ylims=(-3., 3.),
     bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm,
     size=(1000, 1000), linecolor="red", color="red")

plot!(x->x, x->2.75, -3., 3., line=2, color="red", linestyle=:dash, legend=nothing)
savefig("vanderpol_case_1.png")

# --------------------------
# Case 2
# --------------------------

plot_2 = plot(x->x, x->4.0, -2.5, 3., line=2, color="red", linestyle=:dash, legend=nothing)

for i in 1:nsplits_x
    plot!(plot_2, sol_2[i], tickfont=font(30, "Times"), guidefontsize=45,
                   xlab=L"x\raisebox{-0.5mm}{\textcolor{white}{.}}",
                   ylab=L"y\raisebox{2mm}{\textcolor{white}{.}}",
                   xtick=[-2., -1., 0., 1., 2., 3.], ytick=[-4., -3., -2., -1., 0., 1., 2., 3., 4.],
                   xlims=(-2.5, 3.), ylims=(-4.5, 4.),
                   bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm,
                   size=(1000, 1000), color="blue", linewidth=0.0, linecolor="blue", alpha=.5)
end

savefig(plot_2, "vanderpol_case_2.png")

# --------------------------
# Cases 1 and 2 overlapped
# --------------------------

plot_all = plot(x->x, x->4.0, -2.5, 2.5, line=2, color="red", linestyle=:dash, legend=nothing)

for i in 1:nsplits_x
    plot!(plot_all, sol_2[i], tickfont=font(30, "Times"), guidefontsize=45,
                   xlab=L"x\raisebox{-0.5mm}{\textcolor{white}{.}}",
                   ylab=L"y\raisebox{1mm}{\textcolor{white}{.}}",
                   xtick=[-2., -1., 0., 1., 2.], ytick=[-4., -3., -2., -1., 0., 1., 2., 3., 4.],
                   xlims=(-2.5, 2.5), ylims=(-4.5, 4.),
                   bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm,
                   size=(1000, 1000), color="blue", linewidth=0., linecolor="blue")
end

plot!(plot_all, sol_1, tickfont=font(30, "Times"), guidefontsize=45,
                       xlab=L"x\raisebox{-0.5mm}{\textcolor{white}{.}}",
                       ylab=L"y\raisebox{1mm}{\textcolor{white}{.}}",
                       xtick=[-2., -1., 0., 1., 2.], ytick=[-4., -3., -2., -1., 0., 1., 2., 3., 4.],
                       xlims=(-2.5, 2.5), ylims=(-4.5, 4.),
                       bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm,
                       size=(1000, 1000), color="red", linewidth=0., linecolor="red")

plot!(plot_all, x->x, x->2.75, -2.5, 2.5, line=2, color="red", linestyle=:dash, legend=nothing)

savefig(plot_all, "vanderpol_case_all.png")
=#
