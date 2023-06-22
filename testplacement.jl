using Test

include("placement.jl")


@testset "collision function" begin
    
    r = Dict(1 => (Pos(0, 1), Dim(5, 1)), 2 => (Pos(7, 0), Dim(2, 2)))



    # NOK

    @test collision(Pos(0, 0), Dim(2, 3), r)
    @test collision(Pos(0, 1), Dim(1, 5), r)

    @test !collision(Pos(1, 0), Dim(1, 1), r)
    @test collision(Pos(1, 0), Dim(1, 2), r)

    # OK
    @test !collision(Pos(0, 0), Dim(5, 1), r)
    @test !collision(Pos(9, 0), Dim(1, 5), r)

    @test collision(Pos(0, 0), Dim(5, 1), r) || !collision(Pos(0, 0), Dim(0, 1), r)
    @test !collision(Pos(7, 2), Dim(1, 1), r)

    r = Dict( 5 => (Pos(4, 0), Dim(1, 1)), 2 => (Pos(0, 1), Dim(1, 1)), 3 => (Pos(1, 0), Dim(3, 3)), 1 => (Pos(0, 0), Dim(1, 1)))
    @test collision(Pos(0, 2), Dim(2, 1), r)

end

function testoutofbound(r, W)
    @testset "out of bound" begin
        for k in keys(r)
            @test !outofbound(r[k][1], r[k][2], W)
        end
    end
end

@testset "place function" begin
    W = 3
    
    r = place([Dim(1, 1), Dim(1, 1), Dim(2, 1), Dim(1, 1)], W)
    @testset "no collision 1" begin
        for k in keys(r)
            @test !collision(r[k][1], r[k][2], filter(p -> p[1] != k, r))
        end
    end

    testoutofbound(r, W)

    r = place([Dim(1, 1), Dim(1, 1), Dim(2, 1), Dim(1, 1), Dim(7, 1)], W)
    @testset "no collision 2" begin
        for k in keys(r)
            @test !collision(r[k][1], r[k][2], filter(p -> p[1] != k, r))
        end
    end
    testoutofbound(r, W)

    r = place([Dim(7, 1), Dim(1, 1), Dim(1, 1), Dim(2, 1), Dim(1, 1)], W)
    @testset "no collision 3" begin
        for k in keys(r)
            @test !collision(r[k][1], r[k][2], filter(p -> p[1] != k, r))
        end
    end
    testoutofbound(r, W)

    r = place([Dim(1, 1), Dim(3, 1), Dim(1, 1), Dim(2, 1), Dim(1, 1)], W)
    @testset "no collision 4" begin
        for k in keys(r)
            @test !collision(r[k][1], r[k][2], filter(p -> p[1] != k, r))
        end
    end
    testoutofbound(r, W)

    r = place([Dim(1, 1), Dim(1, 1), Dim(3, 3), Dim(2, 1), Dim(1, 1)], W)
    @testset "no collision 5" begin
        for k in keys(r)
            @test !collision(r[k][1], r[k][2], filter(p -> p[1] != k, r))
        end
    end
    testoutofbound(r, W)

    r = place([Dim(1, 1), Dim(4, 1), Dim(2, 1), Dim(8, 2), Dim(1, 1), Dim(3, 1), Dim(4, 2), Dim(2, 1), Dim(1, 1)], W)
    @testset "no collision 6" begin
        for k in keys(r)
            @test !collision(r[k][1], r[k][2], filter(p -> p[1] != k, r))
        end
    end
    testoutofbound(r, W)

end

@testset "genS3 function" begin

    # genS3(W, L, eps)
    maxW = 1
    maxL = 1
    maxeps = 0.1
    for w in 1:maxW
        for le in 1:maxL
            for e in 0.1:0.1:maxeps
                S = genS3(w, le, e)
                @testset "no collision" begin
                    for (k, s) in enumerate(S)
                        r = dict(i => (s[1], s[2]) for (i, s) in enumerate(S))
                        @test !collision(r[k][1], r[k][2], filter(p -> p[1] != k, r))
                    end
                end
                testoutofbound(r, W)
                
            end
        end
    end
end