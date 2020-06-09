#=============================================
install and precompile required Julia packages
=============================================#

# --- installation ---

# install released packages with their exact versions
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

# install unreleased packages
Pkg.add([
    PackageSpec(url="https://github.com/JuliaReach/ReachabilityAnalysis.jl",
                rev="04241504ada5f0fa73f3188c44200fcf14edbcdc"),
    PackageSpec(url="https://github.com/JuliaReach/LazySets.jl",
                rev="7ac9953de8e17095efb5a71d75d5b5dbf61a5d70"), # union sets branch
    PackageSpec(url="https://github.com/JuliaIntervals/.jl",
                rev="1f1f8374941b6ad63699e1563187d6a18558de95") # non-recursive powers branch
])

# --- precompilation ---

import ReachabilityAnalysis
import BenchmarkTools
import Plots
import GR
import LaTeXStrings
