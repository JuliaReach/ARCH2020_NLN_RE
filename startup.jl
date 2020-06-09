# deactivate plot GUI, which is not available in Docker
ENV["GKSwstype"] = "100"

# instantiate project
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

function main()
    println("Running NLN benchmarks...")

    # Van der Pol benchmark
    println("###\nRunning Van der Pol benchmark\n###")
    include("VanDerPol/vanderpol_benchmark.jl")

    # Laub-Loomis benchmark
    println("###\nRunning Laub-Loomis benchmark\n###")
    include("LaubLoomis/laubloomis_benchmark.jl")

    # Quadrotor benchmark
    println("###\nRunning Quadrotor benchmark\n###")
    include("Quadrotor/quadrotor_benchmark.jl")

    println("Finished running benchmarks.")
    nothing
end

main()
