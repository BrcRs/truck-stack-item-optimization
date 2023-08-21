import Pkg
# Pkg.add(["Test", "HTTP", "JSON"])
Pkg.add(["Test", "Coverage", "AutoHashEquals", "Plots", "Documenter"])
# using OAuth, HTTP, JSON
using Test

@testset "testplacement.jl" begin

    include("testplacement.jl")
end

@testset "testordered_stacks.jl" begin
    include("testordered_stacks.jl")
end

@testset "testplacement_visualizer.jl" begin
    # include("testplacement_visualizer.jl")
    # Needs manual testing
end