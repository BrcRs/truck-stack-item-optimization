using Test

include("../src/ordered_stacks.jl")
include("../src/placement_visualizer.jl")
@testset "order_instance(instance::Dict{T, Stack})" begin
    
    instance = Dict(
        1 => Stack(Pos(1, 5), Dim(9, 1)),
        2 => Stack(Pos(9, 5), Dim(6, 4)),
        3 => Stack(Pos(1, 5), Dim(4, 64)),
        4 => Stack(Pos(4, 1), Dim(4, 9))
    )

    res = order_instance(instance)

    sorted_res = sort(res, by= p -> (get_supplier_order(p[2]), get_supplier_dock_order(p[2]), get_plant_dock_order(p[2])))

    @test res == sorted_res

    instance = Dict(
        k => Stack(Pos(rand() * 20, rand() * 20), Dim(rand() * 20, rand() * 20)) for k in 1:10
    )


    res = order_instance(instance)

    sorted_res = sort(res, by= p -> (get_supplier_order(p[2]), get_supplier_dock_order(p[2]), get_plant_dock_order(p[2])))


    @test res == sorted_res

end

@testset "leq_order(s1, s2)" begin
    # leq_order(s1, s2)
    a = OrderedStack(Pos(1, 1), Dim(1, 1), 1, 2, 3)
    b = OrderedStack(Pos(1, 1), Dim(1, 1), 1, 2, 3)

    @test leq_order(a, b)

    a = OrderedStack(Pos(1, 1), Dim(1, 1), 1, 2, 3)
    b = OrderedStack(Pos(1, 1), Dim(1, 1), 1, 2, 2)

    @test !leq_order(a, b)

    a = OrderedStack(Pos(1, 1), Dim(1, 1), 1, 2, 2)
    b = OrderedStack(Pos(1, 1), Dim(1, 1), 1, 2, 3)

    @test leq_order(a, b)

    a = OrderedStack(Pos(1, 1), Dim(1, 1), 1, 1, 3)
    b = OrderedStack(Pos(1, 1), Dim(1, 1), 1, 2, 3)

    @test leq_order(a, b)

    a = OrderedStack(Pos(1, 1), Dim(1, 1), 1, 2, 3)
    b = OrderedStack(Pos(1, 1), Dim(1, 1), 1, 1, 3)

    @test !leq_order(a, b)

    a = OrderedStack(Pos(1, 1), Dim(1, 1), 2, 2, 3)
    b = OrderedStack(Pos(1, 1), Dim(1, 1), 1, 2, 3)

    @test !leq_order(a, b)

    a = OrderedStack(Pos(1, 1), Dim(1, 1), 1, 2, 3)
    b = OrderedStack(Pos(1, 1), Dim(1, 1), 2, 2, 3)

    @test leq_order(a, b)


end

@testset "can_be_placed" begin
    # can_be_placed(solution, o, s::OrderedStack, W, orientation::Symbol; precision=3)
    solution = Dict(
        1 => OrderedStack(Pos(0, 0), Dim(1, 1), 1, 1, 1),
        2 => OrderedStack(Pos(1, 0), Dim(1, 1), 1, 1, 1),
        3 => OrderedStack(Pos(2, 0), Dim(1, 1), 1, 1, 2),
        4 => OrderedStack(Pos(3, 0), Dim(1, 1), 1, 2, 2),
        5 => OrderedStack(Pos(4, 0), Dim(1, 1), 1, 2, 3),
        6 => OrderedStack(Pos(5, 0), Dim(1, 1), 1, 3, 1)
    )

    @test can_be_placed(solution, Pos(4, 1), OrderedStack(Pos(0, 0), Dim(1, 1), 1, 2, 3), 10, :Parallel, verbose=false)
    @test can_be_placed(solution, Pos(4, 1), OrderedStack(Pos(0, 0), Dim(1, 1), 1, 2, 4), 10, :Parallel, verbose=false)
    @test can_be_placed(solution, Pos(4, 1), OrderedStack(Pos(0, 0), Dim(1, 1), 1, 2, 5), 10, :Parallel, verbose=false)
    @test can_be_placed(solution, Pos(4, 1), OrderedStack(Pos(0, 0), Dim(1, 1), 1, 3, 1), 10, :Parallel, verbose=false)

    @test !can_be_placed(solution, Pos(4, 1), OrderedStack(Pos(0, 0), Dim(1, 1), 1, 3, 2), 10, :Parallel, verbose=false)
    @test !can_be_placed(solution, Pos(4, 1), OrderedStack(Pos(0, 0), Dim(1, 1), 1, 3, 3), 10, :Parallel, verbose=false)

    solution = Dict(
        1 => OrderedStack(Pos(5, 0), Dim(1, 1), 1, 1, 1),
        2 => OrderedStack(Pos(5, 0), Dim(1, 1), 1, 1, 2),
        3 => OrderedStack(Pos(5, 0), Dim(1, 1), 1, 1, 3),
        4 => OrderedStack(Pos(5, 0), Dim(1, 1), 1, 1, 4),
        5 => OrderedStack(Pos(5, 0), Dim(1, 1), 1, 1, 5),
        6 => OrderedStack(Pos(5, 0), Dim(1, 1), 1, 1, 6),
    )
    @test !can_be_placed(solution, Pos(2, 0), OrderedStack(Pos(0, 0), Dim(1, 1), 1, 1, 2), 10, :Parallel, verbose=false)

end


