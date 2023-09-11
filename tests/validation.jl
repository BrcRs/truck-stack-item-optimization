import Pkg
# Pkg.add(["Test", "HTTP", "JSON"])
Pkg.add(["Test", "Coverage", "AutoHashEquals", "Plots", "Documenter"])
# using OAuth, HTTP, JSON
using Test

@testset "testplacement.jl" begin

    include("testplacement.jl")
end

@testset "testprojected_pos.jl" begin
    include("testprojected_pos.jl")
end

@testset "testplacement_algorithms.jl" begin
    include("testplacement_algorithms.jl")
end

@testset "testordered_stacks.jl" begin
    include("testordered_stacks.jl")
end

@testset "testplacement_visualizer.jl" begin
    # include("testplacement_visualizer.jl")
    # Needs manual testing
end