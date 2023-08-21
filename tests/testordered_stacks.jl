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

    sorted_res = sort(res, by= p -> (supplier_order(p[2]), supplier_dock_order(p[2]), plant_dock_order(p[2])))

    @test res == sorted_res

    instance = Dict(
        k => Stack(Pos(rand() * 20, rand() * 20), Dim(rand() * 20, rand() * 20)) for k in 1:10
    )


    res = order_instance(instance)

    sorted_res = sort(res, by= p -> (supplier_order(p[2]), supplier_dock_order(p[2]), plant_dock_order(p[2])))


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



@testset "BLtruck loading orders" begin
    instances_solutions = []

    ITER = 1000
    L = 100
    W = 10

    NBCUTS = 20
    NBFUSE = 20
    # most, least = nb_cuts_fuse_calc(5, 20)
    
    # NBCUTS, NBFUSE = most
    
    # NBCUTS = nb_cuts_fuse_avg(4)
    # NBFUSE = 0
    for i in 1:ITER
        rectangles = cutandfuse_generator(L, W, NBCUTS, NBFUSE; precision=3)

        rectangles = order_instance(rectangles)

        instance = shuffle(rectangles)
        solution = nothing

        try
            solution = BLtruck(instance, W; precision=3, loading_order=true)
        catch e
            display(instance)
            throw(e)
        end
        # display(solution)
        foundL = max([get_pos(solution[k]).x + get_dim(solution[k]).le for k in keys(solution)]...)
        push!(instances_solutions, (ins=instance, sol=solution))
    end

    # Test solution is valid
    @testset "BLtruck is valid" begin
        for i in 1:ITER
            instance, solution = instances_solutions[i]
            # check overlapping and out of bounds
            @testset "No overlapping" begin
                for (j, stack) in solution
                    @test !collision(get_pos(stack), get_dim(stack), filter(p -> p[1] != j, solution); precision=3, verbose=false)
                end
            end
            @testset "Not out of bounds" begin
                
                for (j, stack) in solution
                    @test !outofbound(get_pos(stack), get_dim(stack), W; precision=3)
                end
            end
        end 
    end

    ### Not relevant with this constraint
    # # Test optimality is 2-OPT
    # @testset "BLtruck is 2-OPT" begin
    #     for (i, r) in enumerate(ratios)
    #         # @test r <= 2.0 + 0.1
    #         @test r <= 3.0 # TODO why not 2-OPT?


    #         # if r > 2.0 + 0.1
    #         #     display(instances_solutions[i].ins)
    #         #     display(instances_solutions[i].sol)
    #         #     display(length(instances_solutions[i].sol))
    #         #     # println(r)
    #         # end
    #     end
    # end

    # Test loading orders are satisfied
    @testset "loading orders" begin
        for i in 1:ITER
            # make sure the loading order is coherent
            instance, solution = instances_solutions[i]

            x_sorted_sol = sort(collect(solution), by=p -> get_pos(p[2]).x)
            ordered_sol = sort(x_sorted_sol, by= p -> (supplier_order(p[2]), supplier_dock_order(p[2]), plant_dock_order(p[2])))

            # @test [p[1] for p in x_sorted_sol] == [p[1] for p in ordered_sol]
            for p in solution
                @test can_be_placed(filter(x -> x[1] != p[1], solution), get_pos(p[2]), p[2], W, :Parallel)
                if !can_be_placed(filter(x -> x[1] != p[1], solution), get_pos(p[2]), p[2], W, :Parallel)
                    can_be_placed(filter(x -> x[1] != p[1], solution), get_pos(p[2]), p[2], W, :Parallel; verbose=true)
                    println("DEBUG ==========~~~~~~~~------")
                    println(p) 
                    # 11, OrderedStack(Stack(Pos(86.84943832964126, 5.190958474643771), Dim(43.84212073415997, 1.76161657584896)), 1, 5, 1)
                    # 9 => OrderedStack(Stack(Pos(137.609, 2.92317), Dim(13.723, 2.26779)), 1, 5, 1)
                    # 19 => OrderedStack(Stack(Pos(137.609, 0), Dim(43.8421, 2.92317)), 1, 4, 1)
                    display(solution)
                    plot_placement(W, L, solution; orthonormal=true)
                    error()
                end
            end
        end
    end
end

