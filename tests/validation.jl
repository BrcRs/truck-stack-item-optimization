import Pkg
# Pkg.add(["Test", "HTTP", "JSON"])
Pkg.add(["Test", "Coverage"])
# using OAuth, HTTP, JSON
using Test

@testset "testplacement.jl" begin

    include("testplacement.jl")
end

@test 5 - 2 == 3