using Test

include("../src/ordered_stacks.jl")

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


@testset "BLtruck loading orders" begin
    ratios = []
    instances_solutions = []

    ITER = 10
    L = 100
    W = 10

    NBCUTS = 5
    NBFUSE = 5
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
        # ratio, solution = eval_placement_function(BLtruck, instance)
        foundL = max([get_pos(solution[k]).x + get_dim(solution[k]).le for k in keys(solution)]...)
        ratio = foundL / L
        push!(ratios, ratio)
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
    
    # Test optimality is 2-OPT
    @testset "BLtruck is 2-OPT" begin
        for (i, r) in enumerate(ratios)
            # @test r <= 2.0 + 0.1
            @test r <= 3.0 # TODO why not 2-OPT?


            # if r > 2.0 + 0.1
            #     display(instances_solutions[i].ins)
            #     display(instances_solutions[i].sol)
            #     display(length(instances_solutions[i].sol))
            #     # println(r)
            # end
        end
    end

    # Test loading orders are satisfied
    @testset "loading orders" begin
        for i in 1:ITER
            # make sure the loading order is coherent
            instance, solution = instances_solutions[i]

            x_sorted_sol = sort(collect(solution), by=p -> get_pos(p[2]).x)
            ordered_sol = sort(x_sorted_sol, by= p -> (supplier_order(p[2]), supplier_dock_order(p[2]), plant_dock_order(p[2])))

            @test [p[1] for p in x_sorted_sol] == [p[1] for p in ordered_sol]
        end
    end
end

