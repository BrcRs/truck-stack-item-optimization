import Pkg
# Pkg.add(["Test", "HTTP", "JSON"])
Pkg.add(["Test"])
# using OAuth, HTTP, JSON
using Test

@testset "Testing tests" begin
    @test 1 < 2
    @test 2 < 3
    @test 2 + 2 == 4
end

@test 5 - 2 == 3