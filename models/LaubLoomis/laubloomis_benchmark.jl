using BenchmarkTools, Plots, Plots.PlotMeasures, LaTeXStrings
using BenchmarkTools: minimum, median

SUITE = BenchmarkGroup()
model = "LaubLoomis"
cases = ["W=0.01", "W=0.05", "W=0.1"]
SUITE[model] = BenchmarkGroup()

include("laubloomis.jl")
validation = []

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
push!(validation, property ? 1 : 0)

# width of final box
final_width = ρ(e4, sol_1z[end]) + ρ(-e4, sol_1z[end])
println("width of final box, case $(cases[1]) : $final_width")

# benchmark
SUITE["LaubLoomis"][cases[1]] = @benchmarkable solve($prob, T=$T, alg=$alg)

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
final_width = ρ(e4, sol_2z[end]) + ρ(-e4, sol_2z[end])
println("width of final box, case $(cases[2]): $final_width")

# benchmark
SUITE["LaubLoomis"][cases[2]] = @benchmarkable solve($prob, T=$T, alg=$alg)

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
push!(validation, property ? 1 : 0)

# width of final box
final_width = ρ(e4, sol_3z[end]) + ρ(-e4, sol_3z[end])
println("width of final box, case W = $(cases[3]) : $final_width")

# benchmark
SUITE["LaubLoomis"][cases[3]] = @benchmarkable solve($prob, T=$T, alg=$alg)

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

open("runtimes.csv", "w") do io
   print(io, "JuliaReach\n")
   for (i, c) in enumerate(casos)
        print(io, "$c; $(validation[i]); $(runtimes[c])\n")
   end
   print(io, "End of JuliaReach")
end

#=
# ==============================================================================
# Execute benchmarks and save benchmark results
# ==============================================================================

plot(sol_1,
     tickfont=font(30, "Times"), guidefontsize=45,
     xlab=L"t\raisebox{2.0mm}{\textcolor{white}{.}}",
     ylab=L"x_{4}\raisebox{1.2mm}{\textcolor{white}{.}}",
     xtick=[0., 2., 4., 6., 8., 10., 12., 14., 16., 18., 20.],
     ytick=[2, 2.5, 3, 3.5, 4, 4.5],
     xlims=(0., 20.), ylims=(1.5, 4.5),
     bottom_margin=6mm, left_margin=8mm, right_margin=4mm, top_margin=3mm,
     size=(1000, 1000), linecolor="blue")

plot!(x->x, x->4.5, 0., 20., line=2, color="red", linestyle=:dash, legend=nothing)
savefig("laubloomis_case_1.png")

plot(sol_2,
     tickfont=font(30, "Times"), guidefontsize=45,
     xlab=L"t\raisebox{2.0mm}{\textcolor{white}{.}}",
     ylab=L"x_{4}\raisebox{1.2mm}{\textcolor{white}{.}}",
     xtick=[0., 2., 4., 6., 8., 10., 12., 14., 16., 18., 20.],
     ytick=[2, 2.5, 3, 3.5, 4, 4.5, 5.0],
     xlims=(0., 20.), ylims=(1.5, 5.0),
     bottom_margin=6mm, left_margin=8mm, right_margin=4mm, top_margin=3mm,
     size=(1000, 1000), linecolor="blue")

plot!(x->x, x->5.0, 0., 20., line=2, color="red", linestyle=:dash, legend=nothing)
savefig("laubloomis_case_2.png")

plot(sol_3,
     tickfont=font(30, "Times"), guidefontsize=45,
     xlab=L"t\raisebox{2.0mm}{\textcolor{white}{.}}",
     ylab=L"x_{4}\raisebox{1.2mm}{\textcolor{white}{.}}",
     xtick=[0., 2., 4., 6., 8., 10., 12., 14., 16., 18., 20.],
     ytick=[2, 2.5, 3, 3.5, 4, 4.5, 5.0],
     xlims=(0., 20.), ylims=(1.5, 5.0),
     bottom_margin=6mm, left_margin=8mm, right_margin=4mm, top_margin=3mm,
     size=(1000, 1000), linecolor="blue")

plot!(x->x, x->5.0, 0., 20., line=2, color="red", linestyle=:dash, legend=nothing)
savefig("laubloomis_case_3.png")

plot(sol_case_1, color="red")

plot!(sol_case_2, alpha=0.6, color="green")

plot!(sol_case_3, alpha=0.2,
     tickfont=font(30, "Times"), guidefontsize=45,
     xlab=L"t\raisebox{2.0mm}{\textcolor{white}{.}}",
     ylab=L"x_{4}\raisebox{1.2mm}{\textcolor{white}{.}}",
     xtick=[0., 2., 4., 6., 8., 10., 12., 14., 16., 18., 20.],
     ytick=[2, 2.5, 3, 3.5, 4, 4.5, 5.0],
     xlims=(0., 20.), ylims=(1.5, 5.0),
     bottom_margin=6mm, left_margin=8mm, right_margin=4mm, top_margin=3mm,
     size=(1000, 1000), color="blue")

plot!(x->x, x->5.0, 0., 20., line=2, color="red", linestyle=:dash, legend=nothing)
plot!(x->x, x->4.5, 0., 20., line=2, color="red", linestyle=:dash, legend=nothing)
savefig("laubloomis_case_all.png")

=#
