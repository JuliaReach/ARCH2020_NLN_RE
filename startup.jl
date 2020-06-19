# deactivate plot GUI, which is not available in Docker
ENV["GKSwstype"] = "100"

# instantiate project
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

function main()
    global io = open("runtimes.csv", "w")

    println("Running NLN benchmarks...")

    # Production-destruction benchmark
    println("###\nRunning production-destruction model benchmark\n###")
    include("models/ProdDestruct/prod_destruct_benchmark.jl")

    # Coupled Van der Pol benchmark
    println("###\nRunning Van der Pol benchmark\n###")
    include("models/VanDerPol/vanderpol_benchmark.jl")

    # Laub-Loomis benchmark
    println("###\nRunning Laub-Loomis benchmark\n###")
    include("models/LaubLoomis/laubloomis_benchmark.jl")

    # Quadrotor benchmark
    println("###\nRunning Quadrotor benchmark\n###")
    include("models/Quadrotor/quadrotor_benchmark.jl")


    print(io, "\n")
    println("Finished running benchmarks.")
    close(io)
    nothing
end

main()
