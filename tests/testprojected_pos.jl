
include("../src/projected_pos.jl")

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

    @test to_add == [IntersectionPos(5, 1), IntersectionPos(5, 2), IntersectionPos(5, 3)]


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

    @test to_add == [IntersectionPos(1, 5), IntersectionPos(2, 5), IntersectionPos(3, 5)]




    to_add = []

    c = ProjectedPos(Pos(222.5523, 97.3241), Pos(320.1685, 97.3241), :Horizontal)

    allprojected = [
        ProjectedPos(Pos(293.8116, 96.6831), Pos(293.8116, 97.46255), :Vertical),
    ]

    upd_intersection!(to_add, c, allprojected)

    @test to_add == [IntersectionPos(293.8116, 97.3241)]


end

