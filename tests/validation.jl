import Pkg
# Pkg.add(["Test", "HTTP", "JSON"])
Pkg.add(["Test", "Coverage", "AutoHashEquals", "Plots"])
# using OAuth, HTTP, JSON
using Test

@testset "testplacement.jl" begin

    include("testplacement.jl")
end

@testset "testplacement_visualizer.jl" begin
    # include("testplacement_visualizer.jl")
    # Needs manual testing
end