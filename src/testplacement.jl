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

@testset "overlapY" begin
    apos, adim = Pos(0, 0), Dim(1, 1)
    bpos, bdim = Pos(1, 0), Dim(1, 1)
    cpos, cdim = Pos(0, 1), Dim(1, 1)
    dpos, ddim = Pos(1, 1), Dim(1, 1)
    epos, edim = Pos(0.5, 0), Dim(1, 1)
    fpos, fdim = Pos(0, 0.5), Dim(1, 1)
    gpos, gdim = Pos(0.5, 0.5), Dim(1, 1)
    
    @test overlapY(apos, adim, bpos, bdim)
    @test !overlapY(apos, adim, cpos, cdim)
    @test !overlapY(apos, adim, dpos, ddim)
    @test overlapY(apos, adim, epos, edim)
    @test overlapY(apos, adim, fpos, fdim)
    @test overlapY(apos, adim, gpos, gdim)


    @test overlapY(Pos(0, 0), Dim(2, 2), Pos(2, 0), Dim(1, 5))
    @test overlapY(Pos(2, 0), Dim(1, 5), Pos(0, 0), Dim(2, 2))
    # @test !overlapY(Pos(0, 0), Dim(2, 2), Pos(0, 2), Dim(5, 1))
end

@testset "overlapX" begin
    apos, adim = Pos(0, 0), Dim(1, 1)
    bpos, bdim = Pos(1, 0), Dim(1, 1)
    cpos, cdim = Pos(0, 1), Dim(1, 1)
    dpos, ddim = Pos(1, 1), Dim(1, 1)
    epos, edim = Pos(0.5, 0), Dim(1, 1)
    fpos, fdim = Pos(0, 0.5), Dim(1, 1)
    gpos, gdim = Pos(0.5, 0.5), Dim(1, 1)
    
    @test !overlapX(apos, adim, bpos, bdim)
    @test overlapX(apos, adim, cpos, cdim)
    @test !overlapX(apos, adim, dpos, ddim)
    @test overlapX(apos, adim, epos, edim)
    @test overlapX(apos, adim, fpos, fdim)
    @test overlapX(apos, adim, gpos, gdim)
    
    
    # @test overlapX(Pos(0, 0), Dim(2, 2), Pos(0, 2), Dim(5, 1))
    @test !overlapX(Pos(0, 0), Dim(2, 2), Pos(2, 0), Dim(1, 5))
    @test !overlapX(Pos(2, 0), Dim(1, 5), Pos(0, 0), Dim(2, 2))
    
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

    aboxes = findboxesabove(Pos(0.0, 3), Dict(
        1 => (Pos(0, 0), Dim(2, 3))))
    @test isempty(aboxes)

    aboxes = findboxesabove(Pos(2, 0), Dict(
        1 => (Pos(0, 0), Dim(2, 3))))
    @test !(1 in aboxes)

    
    aboxes = findboxesabove(Pos(2, 0), Dict(
        2 => (Pos(0, 3), Dim(10, 1)),
        1 => (Pos(0, 0), Dim(2, 3))))
    @test 2 in aboxes && !(1 in aboxes)

    p = Pos(0.260, 0.378)
    d = Dim(0.198, 0.574)
    r = Dict(
      2 => (Pos(0.260, 0), Dim(0.368, 0.378)),
      3 => (Pos(0.272, 0.496), Dim(0.165, 0.206)),
      1 => (Pos(0, 0), Dim(0.260, 0.496))
    )
    aboxes = findboxesabove(p, r)
    @test !(3 in aboxes)
    @test !(2 in aboxes)
    @test !(1 in aboxes)


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
    aboxes = findboxesabove(p, r)
    @test !(1 in aboxes)
    @test !(2 in aboxes)
    @test !(3 in aboxes)
    @test !(4 in aboxes)
    @test !(5 in aboxes)
    @test !(6 in aboxes)
    @test !(7 in aboxes)

    r = Dict(1 => (Pos(0, 1), Dim(5, 1)), 2 => (Pos(7, 0), Dim(2, 2)))
    
    @test 1 in findboxesabove(Pos(0, 0), r)


    p = Pos(0.9173920482282241, 0.5956575275540524)
    d = Dim(0.001, 0.12376169444300854)
    r = Dict(
      5 => (Pos(0.917392, 0), Dim(0.082608, 0.302143)),
      4 => (Pos(0.751503, 0), Dim(0.165889, 0.367109)),
      6 => (Pos(0.917392, 0.302143), Dim(0.0269158, 0.511427)),
      7 => (Pos(0.944308, 0.302143), Dim(0.0556922, 0.169016)),
      2 => (Pos(0.185148, 0), Dim(0.465636, 0.790959)),
      8 => (Pos(0.916207, 0.700376), Dim(0.00118471, 0.0242651)),
      3 => (Pos(0.650784, 0), Dim(0.100719, 0.727985)),
      1 => (Pos(0, 0), Dim(0.185148, 0.93093))
    ) 
    aboxes = findboxesabove(p, r)

    @test isempty(aboxes)


    p = Pos(0, 0)
    r = Dict(
      1 => (Pos(0, 0), Dim(1, 1))
    )
    aboxes = findboxesabove(p, r)
    @test 1 in aboxes


    p = Pos(1, 0)
    r = Dict(
      1 => (Pos(0, 0), Dim(1, 1))
    )
    aboxes = findboxesabove(p, r)
    @test !(1 in aboxes)


    p = Pos(0.5, 0)
    r = Dict(
      1 => (Pos(0, 1), Dim(1, 1))
    )
    aboxes = findboxesabove(p, r)
    @test 1 in aboxes


    p = Pos(0.5, 0)
    r = Dict(
      1 => (Pos(0, 5), Dim(1, 1))
    )
    aboxes = findboxesabove(p, r)
    @test 1 in aboxes

    p = Pos(1, 0)
    r = Dict(
      1 => (Pos(0, 0), Dim(1, 1))
    )
    aboxes = findboxesabove(p, r)
    @test !(1 in aboxes)

    p = Pos(-0.5, 0)
    r = Dict(
      1 => (Pos(0, 0), Dim(1, 1))
    )
    aboxes = findboxesabove(p, r)
    @test !(1 in aboxes)


end

@testset "find boxes right" begin
    
    r = Dict(1 => (Pos(0, 1), Dim(5, 1)), 2 => (Pos(7, 0), Dim(2, 2)))

    boxesright = findboxesright(Pos(0, 0), Dim(1, 1), r)

    @test !(1 in boxesright)
    @test 2 in boxesright

    p = Pos(0.6600787433709162, 0.5791527054086574)
    d = Dim(0.290875650604264, 0.2395909303148459)
    r = Dict(
      2 => (Pos(0.660079, 0), Dim(0.267742, 0.579153)),
      3 => (Pos(0.671143, 0.606602), Dim(0.123719, 0.271841)),
      1 => (Pos(0, 0), Dim(0.660079, 0.606602))
    )
    boxesright = findboxesright(p, d, r)
    @test !(1 in boxesright)
    @test !(2 in boxesright)
    @test 3 in boxesright

    p = Pos(0.8682584442589909, 0.9947360724710459)
    d = Dim(0.10542070771055019, 0.005263927528954104)
    r = Dict(
    #   5 => (Pos(0.988235, 0), Dim(0.0117654, 0.48294)),
    #   4 => (Pos(0.868258, 0), Dim(0.119976, 0.353601)),
    #   6 => (Pos(0.988235, 0.48294), Dim(0.0117654, 0.486306)),
      7 => (Pos(0.964295, 0.895591), Dim(0.0357055, 0.100487)),
    #   2 => (Pos(0.550112, 0), Dim(0.141112, 0.961887)),
    #   3 => (Pos(0.691224, 0), Dim(0.177034, 0.401109)),
    #   1 => (Pos(0, 0), Dim(0.550112, 0.992972))
    )
    @test 7 in findboxesright(p, d, r)


    r = Dict(1 => (Pos(0, 0), Dim(1, 1)))
    p = Pos(0.0, 0.0)
    d = Dim(1, 1)
    boxesright = findboxesright(p, d, r)

    @test 1 in boxesright

    r = Dict(1 => (Pos(0.5, 0), Dim(1, 1)))
    p = Pos(0.0, 0.0)
    d = Dim(1, 1)
    boxesright = findboxesright(p, d, r)

    @test 1 in boxesright

    r = Dict(1 => (Pos(0, 0), Dim(1, 1)))
    p = Pos(1.0, 0.0)
    d = Dim(1, 1)
    boxesright = findboxesright(p, d, r)

    @test !(1 in boxesright)

    r = Dict(1 => (Pos(0, 0), Dim(1, 1)))
    p = Pos(0.5, 0.0)
    d = Dim(1, 1)
    boxesright = findboxesright(p, d, r)

    @test !(1 in boxesright)

    r = Dict(1 => (Pos(0.5, 0), Dim(1, 1)))
    p = Pos(0.0, 0.0)
    d = Dim(1, 1)
    boxesright = findboxesright(p, d, r)

    @test 1 in boxesright

    r = Dict(1 => (Pos(0.5, 0.5), Dim(1, 1)))
    p = Pos(0.0, 0.0)
    d = Dim(1, 1)
    boxesright = findboxesright(p, d, r)

    @test 1 in boxesright

    r = Dict(1 => (Pos(0.5, -0.5), Dim(1, 1)))
    p = Pos(0.0, 0.0)
    d = Dim(1, 1)
    boxesright = findboxesright(p, d, r)

    @test 1 in boxesright

    r = Dict(1 => (Pos(0, -0.5), Dim(1, 1)))
    p = Pos(0.0, 0.0)
    d = Dim(1, 1)
    boxesright = findboxesright(p, d, r)

    @test 1 in boxesright

    r = Dict(1 => (Pos(0, 1), Dim(1, 1)))
    p = Pos(0.0, 0.0)
    d = Dim(1, 1)
    boxesright = findboxesright(p, d, r)

    @test !(1 in boxesright)

    r = Dict(1 => (Pos(1, 0.5), Dim(1, 1)))
    p = Pos(0.0, 0.0)
    d = Dim(1, 1)
    boxesright = findboxesright(p, d, r)

    @test 1 in boxesright

end

@testset "genWidth" begin
    
    @test genWidth(Pos(0, 0), 2, 0.5) <= 2
    @test genWidth(Pos(0, 1.5), 2, 0.5) == 0.5
    @test genWidth(Pos(0.945, 0), 2, 0.5) <= 2
    @test genWidth(Pos(0.6354684, 0.945), 1.5, 0.2) <= 1.5 - 0.945
    # @test genWidth(Pos(0.649, 0), 2, 0.5) == 2
    # @test genWidth(Pos(0.123, 0), 2, 0.5) == 2
    # @test genWidth(Pos(0.7845, 0), 2, 0.5) == 2
    # @test genWidth(Pos(0.951, 0), 2, 0.5) == 2
    # @test genWidth(Pos(0.6498, 0), 2, 0.5) == 2


    @test_throws ArgumentError genWidth(Pos(0, 0), 2, 3)
    @test_throws ArgumentError genWidth(Pos(0, 1.5), 2, 1)
    @test_throws ArgumentError genWidth(Pos(0, 2), 2, 0.2)
    
end

@testset "genLength" begin
    @test genLength(Pos(0, 0), 2, 0.5) <= 2
    @test genLength(Pos(1.5, 0), 2, 0.5) == 0.5
    @test genLength(Pos(0, 0.945), 2, 0.5) <= 2
    @test genLength(Pos(0.945, 0.6354684), 1.5, 0.2) <= 1.5 - 0.945
    # @test genLength(Pos(0.649, 0), 2, 0.5) == 2
    # @test genLength(Pos(0.123, 0), 2, 0.5) == 2
    # @test genLength(Pos(0.7845, 0), 2, 0.5) == 2
    # @test genLength(Pos(0.951, 0), 2, 0.5) == 2
    # @test genLength(Pos(0.6498, 0), 2, 0.5) == 2


    @test_throws ArgumentError genLength(Pos(0, 0), 2, 3)
    @test_throws ArgumentError genLength(Pos(1.5, 0), 2, 1)
    @test_throws ArgumentError genLength(Pos(2, 0), 2, 0.2)end

@testset "collision function" begin
    
    r = Dict(1 => (Pos(0, 1), Dim(5, 1)), 2 => (Pos(7, 0), Dim(2, 2)))
    
    
    
    # NOK
    @test collision(Pos(0, 0), Dim(2, 3), r)
    @test collision(Pos(0, 1), Dim(1, 5), r)

    @test !collision(Pos(1, 0), Dim(1, 1), r)
    @test collision(Pos(1, 0), Dim(1, 2), r)

    # OK
    @test !collision(Pos(0, 0), Dim(5, 1), r)
    @test !collision(Pos(9, 0), Dim(1, 5), r; verbose=true)

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


    p = Pos(0.8682584442589909, 0.9947360724710459)
    d = Dim(0.10542070771055019, 0.005263927528954104)
    r = Dict(
    #   5 => (Pos(0.988235, 0), Dim(0.0117654, 0.48294)),
    #   4 => (Pos(0.868258, 0), Dim(0.119976, 0.353601)),
    #   6 => (Pos(0.988235, 0.48294), Dim(0.0117654, 0.486306)),
      7 => (Pos(0.964295, 0.895591), Dim(0.0357055, 0.100487)),
    #   2 => (Pos(0.550112, 0), Dim(0.141112, 0.961887)),
    #   3 => (Pos(0.691224, 0), Dim(0.177034, 0.401109)),
    #   1 => (Pos(0, 0), Dim(0.550112, 0.992972))
    )
    @test collision(p, d, r)
    
    p = Pos(0.7801779598676905, 0.9872828760942663)
    d = Dim(0.10785733636911919, 0.012717123905733652)
    r = Dict(
        #   5 => (Pos(0.922803, 0.996865), Dim(0.077197, 0.00313544)),
        #   4 => (Pos(0.922803, 0.893771), Dim(0.077197, 0.103093)),
          6 => (Pos(0.877791, 0.988419), Dim(0.0450124, 0.0115809)),
        # 2 => (Pos(0.922803, 0), Dim(0.077197, 0.781983)),
        # 3 => (Pos(0.922803, 0.781983), Dim(0.077197, 0.111789)),
        # 1 => (Pos(0, 0), Dim(0.922803, 0.879307))
        )
    @test collision(p, d, r)


    p = Pos(0.9173920482282241, 0.5956575275540524)
    d = Dim(0.001, 0.12376169444300854)
    r = Dict(
    #   5 => (Pos(0.917392, 0), Dim(0.082608, 0.302143)),
    #   4 => (Pos(0.751503, 0), Dim(0.165889, 0.367109)),
      6 => (Pos(0.917392, 0.302143), Dim(0.0269158, 0.511427)),
    #   7 => (Pos(0.944308, 0.302143), Dim(0.0556922, 0.169016)),
    #   2 => (Pos(0.185148, 0), Dim(0.465636, 0.790959)),
    #   8 => (Pos(0.916207, 0.700376), Dim(0.00118471, 0.0242651)),
    #   3 => (Pos(0.650784, 0), Dim(0.100719, 0.727985)),
    #   1 => (Pos(0, 0), Dim(0.185148, 0.93093))
    ) 
    @test collision(p, d, r)


    # p = Pos(0.9926480570550624, 0.30824860456814623)
    # d = Dim(0.001, 0.4104763075001061)
    # r = Dict(
    # #   5 => (Pos(0.99176, 0.72853), Dim(0.000888102, 0.106989)),
    # #   4 => (Pos(0.992648, 0.97499), Dim(0.00735194, 0.0250099)),
    # #   2 => (Pos(0.861778, 0), Dim(0.13087, 0.308249)),
    #   3 => (Pos(0.992648, 0), Dim(0.00735194, 0.97499)),
    # #   1 => (Pos(0, 0), Dim(0.861778, 0.971887))
    # )

    # @test !collision(p, d, r)

    p = Pos(0.9859539875368014, 0.9765595603311082)
    d = Dim(0.009885218625116808, 0.00767268558160028)
    r = Dict(
    #   5 => (Pos(0.979308, 0), Dim(0.0175611, 0.150042)),
      12 => (Pos(0.995494, 0.841467), Dim(0.0013758, 0.146971)),
    #   8 => (Pos(0.996869, 0.454622), Dim(0.00313051, 0.474395)),
    #   1 => (Pos(0, 0), Dim(0.864806, 0.830905)),
    #   6 => (Pos(0.996869, 0), Dim(0.00313051, 0.290461)),
    #   11 => (Pos(0.995494, 0.800829), Dim(0.0013758, 0.0406388)),
    #   9 => (Pos(0.996869, 0.929016), Dim(0.00313051, 0.034796)),
    #   3 => (Pos(0.947844, 0), Dim(0.0118101, 0.164162)),
    #   7 => (Pos(0.996869, 0.290461), Dim(0.00313051, 0.16416)),
    #   4 => (Pos(0.959654, 0), Dim(0.0196547, 0.172211)),
    #   13 => (Pos(0.98648, 0.984232), Dim(0.00901378, 0.0102169)),
    #   2 => (Pos(0.864806, 0), Dim(0.0830375, 0.382982)),
    #   10 => (Pos(0.996869, 0.963812), Dim(0.00313051, 0.0293756))
    )
    @test collision(p, d, r)

    r = Dict(1 => (Pos(0, 0), Dim(1, 1)))
    @test collision(Pos(0.9, 0.9), Dim(1/(10^3), 1/(10^3)), r; verbose=true)

end

@testset "illegalpos" begin
    r = Dict(1 => (Pos(0.0, 0.0), Dim(1, 1)))
 
    @test illegalpos(Pos(0, 0), r)
    @test illegalpos(Pos(0, 0.9), r)
    @test illegalpos(Pos(0.9, 0), r)
    @test illegalpos(Pos(0.9, 0.9), r)

    @test !illegalpos(Pos(0, 1), r)
    @test !illegalpos(Pos(1, 0), r)
    @test !illegalpos(Pos(1, 1), r)
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
    maxW = 5
    maxL = 5
    maxeps = 1
    epsstep = 0.1

    my_rs = Dict()

    ratios = Array{Union{Int8, Missing}}(missing, maxW, maxL, convert(Int64, maxeps/epsstep))
    numbers = Array{Union{Int32, Missing}}(missing, maxW, maxL, convert(Int64, maxeps/epsstep))
    # genS3(W, L, eps)
    for w in 1:maxW
        println("w:$w")
        for le in 1:maxL
            println("\tle:$le")
            for (ei, e) in enumerate(0.1:epsstep:maxeps)
                println("\t\te:$e")
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
                ratios[w, le, ei] = convert(Int8, round(volume/(w*le) * 100))
                numbers[w, le, ei] = length(r)
                if isempty(my_rs) && length(r) <= 6 && volume/(w*le) < 60
                    my_rs = r
                end
                testoutofbound(r, w)
                # Filling test
                # println(volume/(w*le))
                # @test volume/(w*le) > 0.9                
            end
        end
    end
    display(ratios)
    display(numbers)
    display(my_rs)
    # figure out why this has a ratio of less than 60%
#     Dict{Any, Any} with 5 entries:
#   5 => (Pos(0.809869, 0.94915), Dim(0.190131, 0.0508503))
#   4 => (Pos(0.809869, 0.734121), Dim(0.190131, 0.215029))
#   2 => (Pos(0.809869, 0), Dim(0.190131, 0.311604))
#   3 => (Pos(0.809869, 0.311604), Dim(0.190131, 0.422517))
#   1 => (Pos(0, 0), Dim(0.809869, 0.813025))
end