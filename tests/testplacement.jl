using Test

include("../src/placement.jl")


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
        1 => Stack(Pos(0, 0), Dim(2, 3))))
    @test isempty(aboxes)

    aboxes = findboxesabove(Pos(2, 0), Dict(
        1 => Stack(Pos(0, 0), Dim(2, 3))))
    @test !(1 in aboxes)

    
    aboxes = findboxesabove(Pos(2, 0), Dict(
        2 => Stack(Pos(0, 3), Dim(10, 1)),
        1 => Stack(Pos(0, 0), Dim(2, 3))))
    @test 2 in aboxes && !(1 in aboxes)

    p = Pos(0.260, 0.378)
    d = Dim(0.198, 0.574)
    r = Dict(
      2 => Stack(Pos(0.260, 0), Dim(0.368, 0.378)),
      3 => Stack(Pos(0.272, 0.496), Dim(0.165, 0.206)),
      1 => Stack(Pos(0, 0), Dim(0.260, 0.496))
    )
    aboxes = findboxesabove(p, r)
    @test !(3 in aboxes)
    @test !(2 in aboxes)
    @test !(1 in aboxes)


    p = Pos(0.679720634248212, 0.5696261712270595)
    d = Dim(0.2343801985845854, 0.17673141164813494)
    r = Dict(
        5 => Stack(Pos(0.916311, 0.509838), Dim(0.0836887, 0.115192)),
        4 => Stack(Pos(0.916311, 0), Dim(0.0836887, 0.509838)),
        6 => Stack(Pos(0.916311, 0.62503), Dim(0.0836887, 0.192314)),
        7 => Stack(Pos(0.779155, 0.629148), Dim(0.117465, 0.264841)),
        2 => Stack(Pos(0.322131, 0), Dim(0.113921, 0.521509)),
        3 => Stack(Pos(0.436052, 0), Dim(0.48026, 0.413936)),
        1 => Stack(Pos(0, 0), Dim(0.322131, 0.837564)))
    aboxes = findboxesabove(p, r)
    @test !(1 in aboxes)
    @test !(2 in aboxes)
    @test !(3 in aboxes)
    @test !(4 in aboxes)
    @test !(5 in aboxes)
    @test !(6 in aboxes)
    @test !(7 in aboxes)

    r = Dict(1 => Stack(Pos(0, 1), Dim(5, 1)), 2 => Stack(Pos(7, 0), Dim(2, 2)))
    
    @test 1 in findboxesabove(Pos(0, 0), r)


    p = Pos(0.9173920482282241, 0.5956575275540524)
    d = Dim(0.001, 0.12376169444300854)
    r = Dict(
      5 => Stack(Pos(0.917392, 0), Dim(0.082608, 0.302143)),
      4 => Stack(Pos(0.751503, 0), Dim(0.165889, 0.367109)),
      6 => Stack(Pos(0.917392, 0.302143), Dim(0.0269158, 0.511427)),
      7 => Stack(Pos(0.944308, 0.302143), Dim(0.0556922, 0.169016)),
      2 => Stack(Pos(0.185148, 0), Dim(0.465636, 0.790959)),
      8 => Stack(Pos(0.916207, 0.700376), Dim(0.00118471, 0.0242651)),
      3 => Stack(Pos(0.650784, 0), Dim(0.100719, 0.727985)),
      1 => Stack(Pos(0, 0), Dim(0.185148, 0.93093))
    ) 
    aboxes = findboxesabove(p, r)

    @test isempty(aboxes)


    p = Pos(0, 0)
    r = Dict(
      1 => Stack(Pos(0, 0), Dim(1, 1))
    )
    aboxes = findboxesabove(p, r)
    @test 1 in aboxes


    p = Pos(1, 0)
    r = Dict(
      1 => Stack(Pos(0, 0), Dim(1, 1))
    )
    aboxes = findboxesabove(p, r)
    @test !(1 in aboxes)


    p = Pos(0.5, 0)
    r = Dict(
      1 => Stack(Pos(0, 1), Dim(1, 1))
    )
    aboxes = findboxesabove(p, r)
    @test 1 in aboxes


    p = Pos(0.5, 0)
    r = Dict(
      1 => Stack(Pos(0, 5), Dim(1, 1))
    )
    aboxes = findboxesabove(p, r)
    @test 1 in aboxes

    p = Pos(1, 0)
    r = Dict(
      1 => Stack(Pos(0, 0), Dim(1, 1))
    )
    aboxes = findboxesabove(p, r)
    @test !(1 in aboxes)

    p = Pos(-0.5, 0)
    r = Dict(
      1 => Stack(Pos(0, 0), Dim(1, 1))
    )
    aboxes = findboxesabove(p, r)
    @test !(1 in aboxes)


end

@testset "find boxes right" begin
    
    r = Dict(1 => Stack(Pos(0, 1), Dim(5, 1)), 2 => Stack(Pos(7, 0), Dim(2, 2)))

    boxesright = findboxesright(Pos(0, 0), Dim(1, 1), r)

    @test !(1 in boxesright)
    @test 2 in boxesright

    p = Pos(0.6600787433709162, 0.5791527054086574)
    d = Dim(0.290875650604264, 0.2395909303148459)
    r = Dict(
      2 => Stack(Pos(0.660079, 0), Dim(0.267742, 0.579153)),
      3 => Stack(Pos(0.671143, 0.606602), Dim(0.123719, 0.271841)),
      1 => Stack(Pos(0, 0), Dim(0.660079, 0.606602))
    )
    boxesright = findboxesright(p, d, r)
    @test !(1 in boxesright)
    @test !(2 in boxesright)
    @test 3 in boxesright

    p = Pos(0.8682584442589909, 0.9947360724710459)
    d = Dim(0.10542070771055019, 0.005263927528954104)
    r = Dict(
    #   5 => Stack(Pos(0.988235, 0), Dim(0.0117654, 0.48294)),
    #   4 => Stack(Pos(0.868258, 0), Dim(0.119976, 0.353601)),
    #   6 => Stack(Pos(0.988235, 0.48294), Dim(0.0117654, 0.486306)),
      7 => Stack(Pos(0.964295, 0.895591), Dim(0.0357055, 0.100487)),
    #   2 => Stack(Pos(0.550112, 0), Dim(0.141112, 0.961887)),
    #   3 => Stack(Pos(0.691224, 0), Dim(0.177034, 0.401109)),
    #   1 => Stack(Pos(0, 0), Dim(0.550112, 0.992972))
    )
    @test 7 in findboxesright(p, d, r)


    r = Dict(1 => Stack(Pos(0, 0), Dim(1, 1)))
    p = Pos(0.0, 0.0)
    d = Dim(1, 1)
    boxesright = findboxesright(p, d, r)

    @test 1 in boxesright

    r = Dict(1 => Stack(Pos(0.5, 0), Dim(1, 1)))
    p = Pos(0.0, 0.0)
    d = Dim(1, 1)
    boxesright = findboxesright(p, d, r)

    @test 1 in boxesright

    r = Dict(1 => Stack(Pos(0, 0), Dim(1, 1)))
    p = Pos(1.0, 0.0)
    d = Dim(1, 1)
    boxesright = findboxesright(p, d, r)

    @test !(1 in boxesright)

    r = Dict(1 => Stack(Pos(0, 0), Dim(1, 1)))
    p = Pos(0.5, 0.0)
    d = Dim(1, 1)
    boxesright = findboxesright(p, d, r)

    @test !(1 in boxesright)

    r = Dict(1 => Stack(Pos(0.5, 0), Dim(1, 1)))
    p = Pos(0.0, 0.0)
    d = Dim(1, 1)
    boxesright = findboxesright(p, d, r)

    @test 1 in boxesright

    r = Dict(1 => Stack(Pos(0.5, 0.5), Dim(1, 1)))
    p = Pos(0.0, 0.0)
    d = Dim(1, 1)
    boxesright = findboxesright(p, d, r)

    @test 1 in boxesright

    r = Dict(1 => Stack(Pos(0.5, -0.5), Dim(1, 1)))
    p = Pos(0.0, 0.0)
    d = Dim(1, 1)
    boxesright = findboxesright(p, d, r)

    @test 1 in boxesright

    r = Dict(1 => Stack(Pos(0, -0.5), Dim(1, 1)))
    p = Pos(0.0, 0.0)
    d = Dim(1, 1)
    boxesright = findboxesright(p, d, r)

    @test 1 in boxesright

    r = Dict(1 => Stack(Pos(0, 1), Dim(1, 1)))
    p = Pos(0.0, 0.0)
    d = Dim(1, 1)
    boxesright = findboxesright(p, d, r)

    @test !(1 in boxesright)

    r = Dict(1 => Stack(Pos(1, 0.5), Dim(1, 1)))
    p = Pos(0.0, 0.0)
    d = Dim(1, 1)
    boxesright = findboxesright(p, d, r)

    @test 1 in boxesright

end

@testset "find boxes left" begin
    
    r = Dict(1 => Stack(Pos(0, 0), Dim(1, 1)))
    p = Pos(0.0, 0.0)
    d = Dim(1, 1)
    boxesleft = findboxesleft(p, d, r)

    @test 1 in boxesleft

    r = Dict(1 => Stack(Pos(0.5, 0), Dim(1, 1)))
    p = Pos(0.0, 0.0)
    d = Dim(1, 1)
    boxesleft = findboxesleft(p, d, r)

    @test !(1 in boxesleft)

    r = Dict(1 => Stack(Pos(0, 0), Dim(1, 1)))
    p = Pos(1.0, 0.0)
    d = Dim(1, 1)
    boxesleft = findboxesleft(p, d, r)

    @test 1 in boxesleft

    r = Dict(1 => Stack(Pos(0, 0), Dim(1, 1)))
    p = Pos(0.5, 0.0)
    d = Dim(1, 1)
    boxesleft = findboxesleft(p, d, r)

    @test 1 in boxesleft

    r = Dict(1 => Stack(Pos(0.5, 0), Dim(1, 1)))
    p = Pos(0.0, 0.0)
    d = Dim(1, 1)
    boxesleft = findboxesleft(p, d, r)

    @test !(1 in boxesleft)

    r = Dict(1 => Stack(Pos(0.5, 0.5), Dim(1, 1)))
    p = Pos(0.0, 0.0)
    d = Dim(1, 1)
    boxesleft = findboxesleft(p, d, r)

    @test !(1 in boxesleft)

    r = Dict(1 => Stack(Pos(0.5, -0.5), Dim(1, 1)))
    p = Pos(0.0, 0.0)
    d = Dim(1, 1)
    boxesleft = findboxesleft(p, d, r)

    @test !(1 in boxesleft)

    r = Dict(1 => Stack(Pos(0, -0.5), Dim(1, 1)))
    p = Pos(0.0, 0.0)
    d = Dim(1, 1)
    boxesleft = findboxesleft(p, d, r)

    @test 1 in boxesleft

    r = Dict(1 => Stack(Pos(0, 1), Dim(1, 1)))
    p = Pos(0.0, 0.0)
    d = Dim(1, 1)
    boxesleft = findboxesleft(p, d, r)

    @test !(1 in boxesleft)

    r = Dict(1 => Stack(Pos(1, 0.5), Dim(1, 1)))
    p = Pos(0.0, 0.0)
    d = Dim(1, 1)
    boxesleft = findboxesleft(p, d, r)

    @test !(1 in boxesleft)


    r = Dict(7 => Stack(Pos(0, 1.79067), Dim(347.168, 98.2086)))
    p = Pos(347.168, 1.79067)
    d = Dim(652.832, 0.01) # this passes
    boxesleft = findboxesleft(p, d, r; verbose=false)

    @test 7 in boxesleft

    r = Dict(4 => Stack(Pos(0, 1.79067), Dim(347.168, 98.2086)))
    p = Pos(347.168, 1.79067)
    d = Dim(652.832, 0.000744365) # but this doesn't
    # with this instance, the precision needs to be higher
    boxesleft = findboxesleft(p, d, r; verbose=false, precision=4) 

    @test 4 in boxesleft

end

@testset "ProjectedPos" begin
    @testset "is_projected" begin
        @test is_projected(ProjectedPos(Pos(5, 0), Pos(5, 5), :Vertical))
        @test !is_projected(Pos(0, 5))
    end

    @testset "set_pos!(::ProjectedPos)" begin
        p = ProjectedPos(Pos(0, 0), Pos(0, 7), :Vertical)
        set_pos!(p, Pos(0, 5))
        @test p == ProjectedPos(Pos(0, 5), Pos(0, 7), :Vertical)

        p = ProjectedPos(Pos(9, 1), Pos(9, 8), :Vertical)
        set_pos!(p, Pos(9, 5.9))
        @test p == ProjectedPos(Pos(9, 5.9), Pos(9, 8), :Vertical)

        p = ProjectedPos(Pos(40.0, 100.0), Pos(50.0, 100.0), :Horizontal)
        set_pos!(p, Pos(25.0, 100.0))
        @test p == ProjectedPos(Pos(25.0, 100.0), Pos(50.0, 100.0), :Horizontal)

    end

    @testset "is_intersected" begin
        p = ProjectedPos(Pos(0, 0), Pos(0, 10), :Vertical)
        @test is_intersected(p, Stack(Pos(0, 0), Dim(1, 1)))


        p = ProjectedPos(Pos(0, 0), Pos(10, 0), :Horizontal)
        @test is_intersected(p, Stack(Pos(0, 0), Dim(1, 1)))

        p = ProjectedPos(Pos(0, 0), Pos(0, 10), :Vertical)
        @test !is_intersected(p, Stack(Pos(1, 0), Dim(1, 1)))

        p = ProjectedPos(Pos(0, 0), Pos(10, 0), :Horizontal)
        @test !is_intersected(p, Stack(Pos(0, 1), Dim(1, 1)))


        p =  ProjectedPos(Pos(0, 0), Pos(0, 10), :Vertical)
        @test is_intersected(p, Stack(Pos(0, 0), Dim(1, 1)))


    end

    @testset "upd!(corners, o::ProjectedPos, s::AbstractStack; precision=3)" begin
        p = ProjectedPos(Pos(2, 0), Pos(2, 10), :Vertical)

        corners = [Pos(0, 0), p, Pos(1, 5)]
        s = Stack(Pos(1, 0), Dim(5, 5))
        upd!(corners, p, s)

        @test corners == [Pos(0, 0), ProjectedPos(Pos(2, 5), Pos(2, 10), :Vertical), Pos(1, 5)]

        p = ProjectedPos(Pos(0, 2), Pos(10, 2), :Horizontal)

        corners = [p, Pos(1, 0), Pos(2, 5)]
        s = Stack(Pos(0, 0), Dim(2, 8))
        upd!(corners, p, s)

        @test corners == [ProjectedPos(Pos(2, 2), Pos(10, 2), :Horizontal), Pos(1, 0), Pos(2, 5)]


        p = ProjectedPos(Pos(0, 2), Pos(10, 2), :Horizontal)

        corners = [p, Pos(1, 0), Pos(2, 5)]
        s = Stack(Pos(0, 0), Dim(11, 8))
        upd!(corners, p, s)

        @test corners == [Pos(1, 0), Pos(2, 5)]


        p = ProjectedPos(Pos(3, 0), Pos(3, 6), :Vertical)

        corners = [Pos(1, 0), p, Pos(4, 5)]
        stacks = [
            Stack(Pos(0, 0), Dim(11, 1)),
            Stack(Pos(0, 1), Dim(10, 0.3)),
            Stack(Pos(2, 1.3), Dim(8, 0.4)),
            Stack(Pos(1.5, 3), Dim(11, 0.2)),
            ]
        for s in stacks
            upd!(corners, p, s)
        end

        @test corners == [Pos(1, 0), ProjectedPos(Pos(3, 3.2), Pos(3, 6), :Vertical), Pos(4, 5)]



        p = ProjectedPos(Pos(3, 0), Pos(3, 6), :Vertical)

        corners = [Pos(1, 0), p, Pos(4, 5)]
        stacks = [
            Stack(Pos(1.5, 3), Dim(11, 0.2)),
            Stack(Pos(2, 1.3), Dim(8, 0.4)),
            Stack(Pos(0, 1), Dim(10, 0.3)),
            Stack(Pos(0, 0), Dim(11, 1)),
            ]
        """
        4                             :
                                      :
                                      :
                                      :
                       1--------------?------------------------------------...
        3              1--------------:------------------------------------...
                                      :
                                      :
                                      :
                                      :
        2                             :
                                      :
                                      :
                                      :
                                      :
        1                             :
                                      :
                                      :
                                      :
                                      :
        0.5                           :
                                      :
                                      :
                                      :
                                      :
        0         1    .    2         X         4
        """



        for s in stacks
            if is_intersected(p, s; verbose=true)
                upd!(corners, p, s)
            end
        end

        @test corners == [Pos(1, 0), ProjectedPos(Pos(3, 3.2), Pos(3, 6), :Vertical), Pos(4, 5)]


        p = ProjectedPos(Pos(5, 0), Pos(5, 5), :Vertical)

        corners = [p, Pos(7, 10)]
        stacks = [
            Stack(Pos(0, 5), Dim(5, 5)),
            Stack(Pos(5, 0), Dim(5, 2)),

            ]
        for s in stacks
            if is_intersected(p, s; verbose=false)
                upd!(corners, p, s; verbose=true)
            end
        end

        @test corners == [ProjectedPos(Pos(5, 2), Pos(5, 5), :Vertical), Pos(7, 10)]



    end

    @testset "intersections" begin
        # upd_intersection!(to_add, c::ProjectedPos, allprojected::Vector{ProjectedPos})

        to_add = []

        c = ProjectedPos(Pos(5, 0), Pos(5, 5), :Vertical)

        allprojected = [
            ProjectedPos(Pos(0, 1), Pos(7, 1), :Horizontal),
            ProjectedPos(Pos(0, 7), Pos(7, 7), :Horizontal),
            ProjectedPos(Pos(0, 2), Pos(7, 2), :Horizontal),
            ProjectedPos(Pos(0, 3), Pos(7, 3), :Horizontal),
            ProjectedPos(Pos(7, 0), Pos(7, 3), :Vertical),
        ]

        upd_intersection!(to_add, c, allprojected; verbose=false)

        @test to_add == [Pos(5, 1), Pos(5, 2), Pos(5, 3)]


        to_add = []

        c = ProjectedPos(Pos(0, 5), Pos(5, 5), :Horizontal)

        allprojected = [
            ProjectedPos(Pos(1, 0), Pos(1, 7), :Vertical),
            ProjectedPos(Pos(7, 0), Pos(7, 7), :Vertical),
            ProjectedPos(Pos(2, 0), Pos(2, 7), :Vertical),
            ProjectedPos(Pos(3, 0), Pos(3, 7), :Vertical),
            ProjectedPos(Pos(0, 7), Pos(3, 7), :Horizontal),
        ]

        upd_intersection!(to_add, c, allprojected)

        @test to_add == [Pos(1, 5), Pos(2, 5), Pos(3, 5)]




        to_add = []

        c = ProjectedPos(Pos(222.5523, 97.3241), Pos(320.1685, 97.3241), :Horizontal)

        allprojected = [
            ProjectedPos(Pos(293.8116, 96.6831), Pos(293.8116, 97.46255), :Vertical),
        ]

        upd_intersection!(to_add, c, allprojected)

        @test to_add == [Pos(293.8116, 97.3241)]


    end

end

@testset "totheleft" begin
    """
    +---+
    | 1 |
    +---p
    """
    solution = Dict(1 => Stack(Pos(0, 0), Dim(1, 1)))
    p = Pos(1.0, 0.0)
    @test totheleft(p, solution; precision=3) == Pos(1.0, 0.0)

    """
    +---+   
    | 1 |   
    +---+   p
    """
    solution = Dict(1 => Stack(Pos(0, 0), Dim(1, 1)))
    p = Pos(2.0, 0.0)
    @test totheleft(p, solution; precision=3) == ProjectedPos(Pos(1.0, 0.0), p, :Horizontal)

    """
    +---+   
    | 1 |   p
    +---+   
    """
    solution = Dict(1 => Stack(Pos(0, 0), Dim(1, 1)))
    p = Pos(2.0, 0.5)
    @test totheleft(p, solution; precision=3) == ProjectedPos(Pos(1.0, 0.5), p, :Horizontal)

    """
    +---+   p
    | 1 |   
    +---+   
    """
    solution = Dict(1 => Stack(Pos(0, 0), Dim(1, 1)))
    p = Pos(2.0, 1)
    @test totheleft(p, solution; precision=3) == ProjectedPos(Pos(0.0, 1.0), p, :Horizontal)

    """
    p    +---+
         | 1 |   
         +---+   
    """
    solution = Dict(1 => Stack(Pos(1, 0), Dim(1, 1)))
    p = Pos(0.0, 1)
    @test totheleft(p, solution; precision=3) == Pos(0.0, 1.0)

    """
         +---+
         | 1 |   
    p    +---+   
    """
    solution = Dict(1 => Stack(Pos(1, 0), Dim(1, 1)))
    p = Pos(0.0, 0.0)
    @test totheleft(p, solution; precision=3) == Pos(0.0, 0.0)


    solution = Dict(1 => Stack(Pos(0, 5.58863), Dim(0.140078, 2.72428)))
    p = Pos(97.2757, 7.53424)
    @test totheleft(p, solution; precision=3) == ProjectedPos(Pos(0.140078, 7.53424), p, :Horizontal)


end


@testset "tothebottom" begin
    """
    +---+
    | 1 |
    +---p
    """
    solution = Dict(1 => Stack(Pos(0, 0), Dim(1, 1)))
    p = Pos(1.0, 0.0)
    @test tothebottom(p, solution; precision=3) == Pos(1.0, 0.0)

    """
    +---+   
    | 1 |   
    +---+   p
    """
    solution = Dict(1 => Stack(Pos(0, 0), Dim(1, 1)))
    p = Pos(2.0, 0.0)
    @test tothebottom(p, solution; precision=3) == Pos(2.0, 0.0)

    """
    +---+   
    | 1 |   p
    +---+   
    """
    solution = Dict(1 => Stack(Pos(0, 0), Dim(1, 1)))
    p = Pos(2.0, 0.5)
    @test tothebottom(p, solution; precision=3) == ProjectedPos(Pos(2.0, 0), p, :Vertical)

    """
    +---+   p
    | 1 |   
    +---+   
    """
    solution = Dict(1 => Stack(Pos(0, 0), Dim(1, 1)))
    p = Pos(2.0, 1)
    @test tothebottom(p, solution; precision=3) == ProjectedPos(Pos(2.0, 0.0), p, :Vertical)

    """
    p    +---+
         | 1 |   
         +---+   
    """
    solution = Dict(1 => Stack(Pos(1, 0), Dim(1, 1)))
    p = Pos(0.0, 1)
    @test tothebottom(p, solution; precision=3) == ProjectedPos(Pos(0.0, 0.0), p, :Vertical)

    """
         +---+
         | 1 |   
    p    +---+   
    """
    solution = Dict(1 => Stack(Pos(1, 0), Dim(1, 1)))
    p = Pos(0.0, 0.0)
    @test tothebottom(p, solution; precision=3) == Pos(0.0, 0.0)


    """
         +---+
         | 1 |   
         +---+p

    """
    solution = Dict(1 => Stack(Pos(1, 1), Dim(1, 1)))
    p = Pos(2.0, 1.0)
    @test tothebottom(p, solution; precision=3) == ProjectedPos(Pos(2.0, 0.0), p, :Vertical)

    """
         p

         +---+
         | 1 |   
         +---+

    """
    solution = Dict(1 => Stack(Pos(1, 1), Dim(1, 1)))
    p = Pos(1.0, 3.0)
    @test tothebottom(p, solution; precision=3) == ProjectedPos(Pos(1.0, 2.0), p, :Vertical)

end


# @testset "genWidth" begin
    
#     @test genWidth(Pos(0, 0), 2, 0.5) <= 2
#     @test genWidth(Pos(0, 1.5), 2, 0.5) == 0.5
#     @test genWidth(Pos(0.945, 0), 2, 0.5) <= 2
#     @test genWidth(Pos(0.6354684, 0.945), 1.5, 0.2) <= 1.5 - 0.945
#     # @test genWidth(Pos(0.649, 0), 2, 0.5) == 2
#     # @test genWidth(Pos(0.123, 0), 2, 0.5) == 2
#     # @test genWidth(Pos(0.7845, 0), 2, 0.5) == 2
#     # @test genWidth(Pos(0.951, 0), 2, 0.5) == 2
#     # @test genWidth(Pos(0.6498, 0), 2, 0.5) == 2


#     @test_throws ArgumentError genWidth(Pos(0, 0), 2, 3)
#     @test_throws ArgumentError genWidth(Pos(0, 1.5), 2, 1)
#     @test_throws ArgumentError genWidth(Pos(0, 2), 2, 0.2)
    
# end

# @testset "genLength" begin
#     @test genLength(Pos(0, 0), 2, 0.5) <= 2
#     @test genLength(Pos(1.5, 0), 2, 0.5) == 0.5
#     @test genLength(Pos(0, 0.945), 2, 0.5) <= 2
#     @test genLength(Pos(0.945, 0.6354684), 1.5, 0.2) <= 1.5 - 0.945
#     # @test genLength(Pos(0.649, 0), 2, 0.5) == 2
#     # @test genLength(Pos(0.123, 0), 2, 0.5) == 2
#     # @test genLength(Pos(0.7845, 0), 2, 0.5) == 2
#     # @test genLength(Pos(0.951, 0), 2, 0.5) == 2
#     # @test genLength(Pos(0.6498, 0), 2, 0.5) == 2


#     @test_throws ArgumentError genLength(Pos(0, 0), 2, 3)
#     @test_throws ArgumentError genLength(Pos(1.5, 0), 2, 1)
#     @test_throws ArgumentError genLength(Pos(2, 0), 2, 0.2)end

@testset "collision function" begin
    
    r = Dict(1 => Stack(Pos(0, 1), Dim(5, 1)),
    2 => Stack(Pos(7, 0), Dim(2, 2)))
    
    """

    +------------------------+         +---------+
    |                        |         |         |
    1------------------------+         |         |
                                       |         |
                                       2---------+

    """
    
    # NOK
    """
    +---------+
    |         |
    +---------|--------------+         +---------+
    |         |              |         |         |
    1---------|--------------+         |         |
    |         |                        |         |
    +---------+                        2---------+

    """

    @test collision(Pos(0, 0), Dim(2, 3), r)

    """
    ......
    |   |
    |   |
    +---|--------------------+         +---------+
    |   |                    |         |         |
    1===+--------------------+         |         |
                                       |         |
                                       2---------+

    """
    
    @test collision(Pos(0, 1), Dim(1, 5), r)
    
    """

    +------------------------+         +---------+
    |                        |         |         |
    1---+---+----------------+         |         |
        |   |                          |         |
        +---+                          2---------+

    """
    
    @test !collision(Pos(1, 0), Dim(1, 1), r)
    """

    +---+---+----------------+         +---------+
    |   |   |                |         |         |
    1---|---|----------------+         |         |
        |   |                          |         |
        +---+                          2---------+

    """
    
    @test collision(Pos(1, 0), Dim(1, 2), r)

    # OK
    @test !collision(Pos(0, 0), Dim(5, 1), r)
    @test !collision(Pos(9, 0), Dim(1, 5), r; verbose=true)

    @test collision(Pos(0, 0), Dim(5, 1), r) || !collision(Pos(0, 0), Dim(0.001, 1), r)
    @test !collision(Pos(7, 2), Dim(1, 1), r)

    @test !collision(Pos(9, 0), Dim(1, 1), r)

    r = Dict( 5 => Stack(Pos(4, 0), Dim(1, 1)), 2 => Stack(Pos(0, 1), Dim(1, 1)), 3 => Stack(Pos(1, 0), Dim(3, 3)), 1 => Stack(Pos(0, 0), Dim(1, 1)))
    @test collision(Pos(0, 2), Dim(2, 1), r)

    p = Pos(0.2604926233999792, 0.37827530439454293)
    d = Dim(0.19831329972752604, 0.5736318825516512)
    r = Dict(
      2 => Stack(Pos(0.260493, 0), Dim(0.367593, 0.378275))
    )
    @test !collision(p, d, r)

    r = Dict(
      3 => Stack(Pos(0.272358, 0.496099), Dim(0.164996, 0.206335))
    )
    @test collision(p, d, r)

    r = Dict(
      1 => Stack(Pos(0, 0), Dim(0.260493, 0.496099))
    )
    @test !collision(p, d, r)


    p = Pos(0.8682584442589909, 0.9947360724710459)
    d = Dim(0.10542070771055019, 0.005263927528954104)
    r = Dict(
    #   5 => Stack(Pos(0.988235, 0), Dim(0.0117654, 0.48294)),
    #   4 => Stack(Pos(0.868258, 0), Dim(0.119976, 0.353601)),
    #   6 => Stack(Pos(0.988235, 0.48294), Dim(0.0117654, 0.486306)),
      7 => Stack(Pos(0.964295, 0.895591), Dim(0.0357055, 0.100487)),
    #   2 => Stack(Pos(0.550112, 0), Dim(0.141112, 0.961887)),
    #   3 => Stack(Pos(0.691224, 0), Dim(0.177034, 0.401109)),
    #   1 => Stack(Pos(0, 0), Dim(0.550112, 0.992972))
    )
    @test collision(p, d, r)
    
    p = Pos(0.7801779598676905, 0.9872828760942663)
    d = Dim(0.10785733636911919, 0.012717123905733652)
    r = Dict(
        #   5 => Stack(Pos(0.922803, 0.996865), Dim(0.077197, 0.00313544)),
        #   4 => Stack(Pos(0.922803, 0.893771), Dim(0.077197, 0.103093)),
          6 => Stack(Pos(0.877791, 0.988419), Dim(0.0450124, 0.0115809)),
        # 2 => Stack(Pos(0.922803, 0), Dim(0.077197, 0.781983)),
        # 3 => Stack(Pos(0.922803, 0.781983), Dim(0.077197, 0.111789)),
        # 1 => Stack(Pos(0, 0), Dim(0.922803, 0.879307))
        )
    @test collision(p, d, r)


    p = Pos(0.9173920482282241, 0.5956575275540524)
    d = Dim(0.001, 0.12376169444300854)
    r = Dict(
    #   5 => Stack(Pos(0.917392, 0), Dim(0.082608, 0.302143)),
    #   4 => Stack(Pos(0.751503, 0), Dim(0.165889, 0.367109)),
      6 => Stack(Pos(0.917392, 0.302143), Dim(0.0269158, 0.511427)),
    #   7 => Stack(Pos(0.944308, 0.302143), Dim(0.0556922, 0.169016)),
    #   2 => Stack(Pos(0.185148, 0), Dim(0.465636, 0.790959)),
    #   8 => Stack(Pos(0.916207, 0.700376), Dim(0.00118471, 0.0242651)),
    #   3 => Stack(Pos(0.650784, 0), Dim(0.100719, 0.727985)),
    #   1 => Stack(Pos(0, 0), Dim(0.185148, 0.93093))
    ) 
    @test collision(p, d, r)


    # p = Pos(0.9926480570550624, 0.30824860456814623)
    # d = Dim(0.001, 0.4104763075001061)
    # r = Dict(
    # #   5 => Stack(Pos(0.99176, 0.72853), Dim(0.000888102, 0.106989)),
    # #   4 => Stack(Pos(0.992648, 0.97499), Dim(0.00735194, 0.0250099)),
    # #   2 => Stack(Pos(0.861778, 0), Dim(0.13087, 0.308249)),
    #   3 => Stack(Pos(0.992648, 0), Dim(0.00735194, 0.97499)),
    # #   1 => Stack(Pos(0, 0), Dim(0.861778, 0.971887))
    # )

    # @test !collision(p, d, r)

    p = Pos(0.9859539875368014, 0.9765595603311082)
    d = Dim(0.009885218625116808, 0.00767268558160028)
    r = Dict(
    #   5 => Stack(Pos(0.979308, 0), Dim(0.0175611, 0.150042)),
      12 => Stack(Pos(0.995494, 0.841467), Dim(0.0013758, 0.146971)),
    #   8 => Stack(Pos(0.996869, 0.454622), Dim(0.00313051, 0.474395)),
    #   1 => Stack(Pos(0, 0), Dim(0.864806, 0.830905)),
    #   6 => Stack(Pos(0.996869, 0), Dim(0.00313051, 0.290461)),
    #   11 => Stack(Pos(0.995494, 0.800829), Dim(0.0013758, 0.0406388)),
    #   9 => Stack(Pos(0.996869, 0.929016), Dim(0.00313051, 0.034796)),
    #   3 => Stack(Pos(0.947844, 0), Dim(0.0118101, 0.164162)),
    #   7 => Stack(Pos(0.996869, 0.290461), Dim(0.00313051, 0.16416)),
    #   4 => Stack(Pos(0.959654, 0), Dim(0.0196547, 0.172211)),
    #   13 => Stack(Pos(0.98648, 0.984232), Dim(0.00901378, 0.0102169)),
    #   2 => Stack(Pos(0.864806, 0), Dim(0.0830375, 0.382982)),
    #   10 => Stack(Pos(0.996869, 0.963812), Dim(0.00313051, 0.0293756))
    )
    @test collision(p, d, r)

    r = Dict(1 => Stack(Pos(0, 0), Dim(1, 1)))
    @test collision(Pos(0.9, 0.9), Dim(1/(10^3), 1/(10^3)), r; verbose=true)


    
    # ```
    #     +---+
    #     | 1 |
    # +---+---+---+
    # | 3 |   | 2 |
    # +---+---+---+
    #     | 4 |
    #     +---+
    # ```
    
    r = Dict(1 => Stack(Pos(1, 2), Dim(1, 1)), 2 => Stack(Pos(2, 1), Dim(1, 1)), 3 => Stack(Pos(0, 1), Dim(1, 1)), 4 => Stack(Pos(1, 0), Dim(1, 1)), 5 => Stack(Pos(1, 1), Dim(1, 1)))
    
    for k in keys(r)    
        @test !collision(r[k].pos, r[k].dim, filter(p -> p[1] != k, r); verbose=true)
    end

end

@testset "illegalpos" begin
    r = Dict(1 => Stack(Pos(0.0, 0.0), Dim(1, 1)))
 
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
            @test !outofbound(r[k].pos, r[k].dim, W)
        end
    end
end

@testset "coveredcorners" begin
    corners = [Pos(0, 0), Pos(5, 5), Pos(5, 0), Pos(0, 5), Pos(1, 1), Pos(0, 6), Pos(6, 0)]
    stack = Stack(Pos(0, 0), Dim(5, 5))
    o = Pos(0, 0)
    shouldbecovered = [Pos(0, 0), Pos(1, 1)]
    covered = coveredcorners(corners, o, stack.dim.le, stack.dim.wi)
    for c in shouldbecovered
        @test c in covered
    end
    for c in filter(x -> !(x in shouldbecovered), corners)
        @test !(c in covered)
    end

    corners = [ Pos(0, 9.980043825131686), 
                Pos(0.6467206477622178, 9.962824093254572), 
                Pos(0.6467206477622178, 9.962824093254572), 
                Pos(2.8535725859835184, 2.441387331367891), 
                Pos(62.623393342666866, 9.239718358416525), 
                Pos(65.73516828201012, 4.490758002528395), 
                Pos(67.98165946491517, 8.754314549106567), 
                Pos(76.47110582244056, 4.490758002528395), 
                Pos(76.51649610654447, 0), 
                Pos(80.45705471972516, 0)]
    pos = Pos(76.51649610654447, 0)
    s = Stack(Pos(9.512621318536885, 6.888225060656747), Dim(3.1117749393432534, 4.629616629616873))
    precision = 3

    covered = coveredcorners(corners, pos, s.dim.le, s.dim.wi; precision=precision, verbose=true)
    @test !(Pos(80.45705471972516, 0) in covered)
    


end

@testset "is_secure" begin

    stack0 = Stack(Pos(0, 0), Dim(1, 1))
    stack1 = Stack(Pos(1.0, 0.0), Dim(1, 1))
    stack2 = Stack(Pos(2.0, 0.0), Dim(1, 1))
    
    solution = Dict(1 => stack0)
    
    @test is_secure(stack1, solution; precision=3)
    @test !is_secure(stack2, solution; precision=3)

    @test is_secure(stack0, Dict(); precision=3)

end

