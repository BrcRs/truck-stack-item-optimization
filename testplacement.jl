using Test

include("placement.jl")


@testset "comparisons with rounding" begin
    @testset "leqtol" begin
        @test leqtol(1, 1, 0)
        @test leqtol(1, 1, 1)
        @test leqtol(1, 1.1, 0)
        @test leqtol(1, 1.1, 1)
        @test !leqtol(1.1, 1, 1)
        @test leqtol(1.111, 1, 0)
        @test !leqtol(1.111, 1, 1)
        @test leqtol(1.111, 2.111, 0)
        @test leqtol(2.111, 2.112, 3)
        @test !leqtol(2.11166, 2.11111, 3)
        @test !leqtol(2.111, 1.111, 0)
        @test leqtol(2.111, 1.6, 0)
        @test !leqtol(2.0, 1.6, 1)
        # @test leqtol(1, 1, 0)
    end
    
    @testset "greatertol" begin
        @test !greatertol(1, 1, 0)
        @test !greatertol(1, 1, 1)
        @test !greatertol(1, 1.1, 0)
        @test !greatertol(1, 1.1, 1)
        @test greatertol(1.1, 1, 1)
        @test !greatertol(1.111, 1, 0)
        @test greatertol(1.111, 1, 1)
        @test !greatertol(1.111, 2.111, 0)
        @test !greatertol(2.111, 2.112, 3)
        @test greatertol(2.11166, 2.11111, 3)
        @test greatertol(2.111, 1.111, 0)
        @test !greatertol(2.111, 1.6, 0)
        @test greatertol(2.0, 1.6, 1)
        # @test greatertol(1, 1, 0)
    end
    
    @testset "lessertol" begin
        @test !lessertol(1, 1, 0)
        @test !lessertol(1, 1, 1)
        @test !lessertol(1, 1.1, 0)
        @test lessertol(1, 1.1, 1)
        @test !lessertol(1.1, 1, 1)
        @test !lessertol(1.111, 1, 0)
        @test !lessertol(1.111, 1, 1)
        @test lessertol(1.111, 2.111, 0)
        @test lessertol(2.111, 2.112, 3)
        @test !lessertol(2.11166, 2.11111, 3)
        @test !lessertol(2.111, 1.111, 0)
        @test !lessertol(2.111, 1.6, 0)
        @test !lessertol(2.0, 1.6, 1)
        # @test leqtol(1, 1, 0)
    end
    
    @testset "geqtol" begin
        @test geqtol(1, 1, 0)
        @test geqtol(1, 1, 1)
        @test geqtol(1, 1.1, 0)
        @test !geqtol(1, 1.1, 1)
        @test geqtol(1.1, 1, 1)
        @test geqtol(1.111, 1, 0)
        @test geqtol(1.111, 1, 1)
        @test !geqtol(1.111, 2.111, 0)
        @test !geqtol(2.111, 2.112, 3)
        @test geqtol(2.11166, 2.11111, 3)
        @test geqtol(2.111, 1.111, 0)
        @test geqtol(2.111, 1.6, 0)
        @test geqtol(2.0, 1.6, 1)
        # @test leqtol(1, 1, 0)
    end
    
    @testset "eqtol" begin
        @test eqtol(1, 1, 0)
        @test eqtol(1, 1, 1)
        @test eqtol(1, 1.1, 0)
        @test !eqtol(1, 1.1, 1)
        @test !eqtol(1.1, 1, 1)
        @test eqtol(1.111, 1, 0)
        @test !eqtol(1.111, 1, 1)
        @test !eqtol(1.111, 2.111, 0)
        @test !eqtol(2.111, 2.112, 3)
        @test !eqtol(2.11166, 2.11111, 3)
        @test !eqtol(2.111, 1.111, 0)
        @test eqtol(2.111, 1.6, 0)
        @test !eqtol(2.0, 1.6, 1)
        # @test leqtol(1, 1, 0)
    end
    
    
    # greatertol(a, b, decimals=3) = geqtol(a, b+(1/decimals), decimals)
    # lessertol(a, b, decimals=3) = leqtol(a, b-(1/decimals), decimals)
    
    # geqtol(a, b, decimals=3) = round(a, digits=decimals) >= round(b, digits=decimals)
    
    # eqtol(a, b, decimals=3) = round(a, digits=decimals) == round(b, digits=decimals)
end

@testset "findboxesabove" begin
    # aboxes = findboxesabove(Pos(0.8461362603452564, 0.8943021672958141), Dim(0.12541640596916293, 0.10284447972870292), Dict(
    #     5 => (Pos(0.940318, 0.898831), Dim(0.0596824, 0.100995)), 
    #     4 => (Pos(1.0, 0), Dim(0.0, 0.371941)), 
    #     2 => (Pos(0.86661, 0), Dim(0.132198, 0.39929)), 
    #     3 => (Pos(0.998808, 0), Dim(0.00119218, 0.19648)), 
    #     1 => (Pos(0, 0), Dim(0.86661, 0.784772))))
    # display(aboxes)
    # @test !isempty(aboxes)#?

    aboxes = findboxesabove(Pos(0.0, 3), Dim(1, 2), Dict(
        1 => (Pos(0, 0), Dim(2, 3))))
    display(aboxes)
    @test isempty(aboxes)

    aboxes = findboxesabove(Pos(2, 0), Dim(1, 2), Dict(
        1 => (Pos(0, 0), Dim(2, 3))))
    display(aboxes)
    @test 1 in aboxes

    
    aboxes = findboxesabove(Pos(2, 0), Dim(1, 1), Dict(
        2 => (Pos(0, 3), Dim(10, 1)),
        1 => (Pos(0, 0), Dim(2, 3))))
    display(aboxes)
    @test 2 in aboxes && 1 in aboxes

    p = Pos(0.2604926233999792, 0.37827530439454293)
    d = Dim(0.19831329972752604, 0.5736318825516512)
    r = Dict(
      2 => (Pos(0.260493, 0), Dim(0.367593, 0.378275)),
      3 => (Pos(0.272358, 0.496099), Dim(0.164996, 0.206335)),
      1 => (Pos(0, 0), Dim(0.260493, 0.496099))
    )
    aboxes = findboxesabove(p, d, r)
    @test 3 in aboxes && !(2 in aboxes) && !(1 in aboxes)


    p = Pos(0.679720634248212, 0.5696261712270595)
    d = Dim(0.2343801985845854, 0.17673141164813494)
    r = Dict(
        5 => (Pos(0.916311, 0.509838), Dim(0.0836887, 0.115192)),
        4 => (Pos(0.916311, 0), Dim(0.0836887, 0.509838)),
        6 => (Pos(0.916311, 0.62503), Dim(0.0836887, 0.192314)),
        7 => (Pos(0.779155, 0.629148), Dim(0.117465, 0.264841)),
        2 => (Pos(0.322131, 0), Dim(0.113921, 0.521509)),
        3 => (Pos(0.436052, 0), Dim(0.48026, 0.413936)),
        1 => (Pos(0, 0), Dim(0.322131, 0.837564)))
    aboxes = findboxesabove(p, d, r)
    @test 7 in aboxes

    aboxes = findboxesabove(p, Dim(0.001, 0.001), r)
    @test 7 in aboxes

    r = Dict(1 => (Pos(0, 1), Dim(5, 1)), 2 => (Pos(7, 0), Dim(2, 2)))
    
    @test 1 in findboxesabove(Pos(0, 0), Dim(2, 3), r)

end

@testset "find boxes right" begin
    
    r = Dict(1 => (Pos(0, 1), Dim(5, 1)), 2 => (Pos(7, 0), Dim(2, 2)))


    @test 1 in findboxesright(Pos(0, 0), Dim(5, 1), r)

    p = Pos(0.6600787433709162, 0.5791527054086574)
    d = Dim(0.290875650604264, 0.2395909303148459)
    r = Dict(
      2 => (Pos(0.660079, 0), Dim(0.267742, 0.579153)),
      3 => (Pos(0.671143, 0.606602), Dim(0.123719, 0.271841)),
      1 => (Pos(0, 0), Dim(0.660079, 0.606602))
    )

    @test 3 in findboxesright(p, d, r)


end

@testset "collision function" begin
    
    r = Dict(1 => (Pos(0, 1), Dim(5, 1)), 2 => (Pos(7, 0), Dim(2, 2)))
    
    
    
    # NOK
    @test collision(Pos(0, 0), Dim(2, 3), r)
    @test collision(Pos(0, 1), Dim(1, 5), r)

    @test !collision(Pos(1, 0), Dim(1, 1), r)
    @test collision(Pos(1, 0), Dim(1, 2), r)

    # OK
    display(findboxesabove(Pos(0, 0), Dim(5, 1), r))
    display(findboxesright(Pos(0, 0), Dim(5, 1), r))
    @test !collision(Pos(0, 0), Dim(5, 1), r)
    @test !collision(Pos(9, 0), Dim(1, 5), r)

    @test collision(Pos(0, 0), Dim(5, 1), r) || !collision(Pos(0, 0), Dim(0.001, 1), r)
    @test !collision(Pos(7, 2), Dim(1, 1), r)

    @test !collision(Pos(9, 0), Dim(1, 1), r)

    r = Dict( 5 => (Pos(4, 0), Dim(1, 1)), 2 => (Pos(0, 1), Dim(1, 1)), 3 => (Pos(1, 0), Dim(3, 3)), 1 => (Pos(0, 0), Dim(1, 1)))
    @test collision(Pos(0, 2), Dim(2, 1), r)

    p = Pos(0.2604926233999792, 0.37827530439454293)
    d = Dim(0.19831329972752604, 0.5736318825516512)
    r = Dict(
      2 => (Pos(0.260493, 0), Dim(0.367593, 0.378275))
    )
    @test !collision(p, d, r)

    r = Dict(
      3 => (Pos(0.272358, 0.496099), Dim(0.164996, 0.206335))
    )
    @test collision(p, d, r)

    r = Dict(
      1 => (Pos(0, 0), Dim(0.260493, 0.496099))
    )
    @test !collision(p, d, r)
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
    maxW = 5
    maxL = 5
    maxeps = 1
    for w in 1:maxW
        # println("w:$w")
        for le in 1:maxL
            # println("\tle:$le")
            for e in 0.1:0.1:maxeps
                # println("\t\te:$e")
                volume = 0.0
                S, r = genS3(w, le, e)
                @testset "no collision" begin
                    for k in keys(r)
                        # r = Dict(
                        #     i => (s.le, s.wi) for (i, s) in enumerate(S))
                        @test !collision(r[k][1], r[k][2], filter(p -> p[1] != k, r))
                        volume += r[k][2].le * r[k][2].wi
                    end
                end
                testoutofbound(r, w)
                # Filling test
                @test volume/(w*le) > 0.9                
            end
        end
    end
end