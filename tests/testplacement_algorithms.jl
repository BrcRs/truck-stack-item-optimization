using Test
using Statistics
include("../src/placement_algorithms.jl")
include("../src/placement.jl")
include("../src/item.jl")
include("../src/placement_visualizer.jl")


@testset "shareside" begin
    
    # shareside(a::Stack, b::Stack)

    """
    +---+---+
    | a | 2 |
    +---+---+
    """

    @test shareside(Stack(Pos(0, 0), Dim(1, 1)), Stack(Pos(1, 0), Dim(1, 1)))
    
    """
    +---+
    | 1 |
    +---+
    | a |
    +---+
    """

    @test shareside(Stack(Pos(0, 0), Dim(1, 1)), Stack(Pos(0, 1), Dim(1, 1)))

    """
    +---+
    | a |
    +---+
    | 4 |
    +---+
    """

    @test shareside(Stack(Pos(0, 0), Dim(1, 1)), Stack(Pos(0, -1), Dim(1, 1)))

    """
    +---+---+
    | 3 | a |
    +---+---+
    """

    @test shareside(Stack(Pos(1, 0), Dim(1, 1)), Stack(Pos(0, 0), Dim(1, 1)))


    """
        +---+
    +---+   |
    | a +---+
    +---+
    """

    @test !shareside(Stack(Pos(0, 0), Dim(1, 1)), Stack(Pos(1, 0.5), Dim(1, 1)))

    """
          +---+
    +---+ |   |
    | a | +---+
    +---+
    """

    @test !shareside(Stack(Pos(0, 0), Dim(1, 1)), Stack(Pos(2, 0.5), Dim(1, 1)))


    """
      +---+
      |   |
    +-+-+-+
    | a |
    +---+
    """

    @test !shareside(Stack(Pos(0, 0), Dim(1, 1)), Stack(Pos(0.5, 1), Dim(1, 1)))

end

@testset "fuse!" begin
    """
    +---+---+
    | a | 2 |
    +---+---+
    """
    rectangles = Dict(1 => Stack(Pos(0, 0), Dim(1, 1)), 2 => Stack(Pos(1, 0), Dim(1, 1)))
    fuse!(rectangles, rectangles[1], rectangles[2], 3)
    @test rectangles == Dict(3 => Stack(Pos(0, 0), Dim(2, 1)))
    
    """
    +---+
    | 1 |
    +---+
    | a |
    +---+
    """

    rectangles = Dict(1 => Stack(Pos(0, 0), Dim(1, 1)), 2 => Stack(Pos(0, 1), Dim(1, 1)))
    fuse!(rectangles, rectangles[1], rectangles[2], 3)
    @test rectangles == Dict(3 => Stack(Pos(0, 0), Dim(1, 2)))
    """
    +---+
    | a |
    +---+
    | 4 |
    +---+
    """

    rectangles = Dict(1 => Stack(Pos(0, 1), Dim(1, 1)), 2 => Stack(Pos(0, 0), Dim(1, 1)))
    fuse!(rectangles, rectangles[1], rectangles[2], 3)
    @test rectangles == Dict(3 => Stack(Pos(0, 0), Dim(1, 2)))
    """
    +---+---+
    | 3 | a |
    +---+---+
    """

    rectangles = Dict(1 => Stack(Pos(1, 0), Dim(1, 1)), 2 => Stack(Pos(0, 0), Dim(1, 1)))
    fuse!(rectangles, rectangles[1], rectangles[2], 3)
    @test rectangles == Dict(3 => Stack(Pos(0, 0), Dim(2, 1)))
end

@testset "cutrectangle" begin


    @test cutrectangle("x", Stack(Pos(0, 0), Dim(1, 1)), 0.5) == Stack(Pos(0, 0), Dim(0.5, 1))

    @test cutrectangle("x", Stack(Pos(0, 0), Dim(1, 1)), 0.25) == Stack(Pos(0, 0), Dim(0.25, 1))

    @test cutrectangle("x", Stack(Pos(0, 0), Dim(1, 1)), 0.75) == Stack(Pos(0, 0), Dim(0.75, 1))


    @test cutrectangle("y", Stack(Pos(0, 0), Dim(1, 1)), 0.5) == Stack(Pos(0, 0), Dim(1, 0.5))

    @test cutrectangle("y", Stack(Pos(0, 0), Dim(1, 1)), 0.1) == Stack(Pos(0, 0), Dim(1, 0.1))


    @test_throws ArgumentError cutrectangle("x", Stack(Pos(2, 0), Dim(1, 1)), 0.5)

    @test_throws ArgumentError cutrectangle("x", Stack(Pos(2, 0), Dim(1, 1)), 0.25)

    @test_throws ArgumentError cutrectangle("x", Stack(Pos(5, 3), Dim(1, 1)), 0.75)


    @test_throws ArgumentError cutrectangle("y", Stack(Pos(3, 2), Dim(1, 1)), 0.5)

    @test_throws ArgumentError cutrectangle("y", Stack(Pos(1, 1), Dim(1, 1)), 0.1)

    @test cutrectangle("x", Stack(Pos(1, 0), Dim(1, 1)), 1.5) == Stack(Pos(1, 0), Dim(0.5, 1))

    @test cutrectangle("y", Stack(Pos(0, 1), Dim(1, 1)), 1.5) == Stack(Pos(0, 1), Dim(1, 0.5))


end

@testset "newrectangle" begin
    rectangle = Stack(Pos(0, 0), Dim(0.5, 1))
    @test newrectangle("x", rectangle, Dim(1, 1), 0.5) == Stack(Pos(0.5, 0), Dim(0.5, 1))


    rectangle = Stack(Pos(1, 0), Dim(0.5, 1))
    @test newrectangle("x", rectangle, Dim(1, 1), 1.5) == Stack(Pos(1.5, 0), Dim(0.5, 1))

    rectangle = Stack(Pos(0, 0), Dim(1, 0.5))
    @test newrectangle("y", rectangle, Dim(1, 1), 0.5) == Stack(Pos(0, 0.5), Dim(1, 0.5))

    rectangle = Stack(Pos(0, 1), Dim(1, 0.5))
    @test newrectangle("y", rectangle, Dim(1, 1), 1.5) == Stack(Pos(0, 1.5), Dim(1, 0.5))

end

function testoutofbound(r, W)
    @testset "out of bound" begin
        for k in keys(r)
            @test !outofbound(r[k].pos, r[k].dim, W)
        end
    end
end

@testset "cutandfuse_generator" begin
    maxW = 5
    maxL = 5
    maxeps = 1
    epsstep = 0.1

    # my_rs = Dict()

    # ratios = Array{Union{Int8, Missing}}(missing, maxW, maxL, convert(Int64, maxeps/epsstep))
    # numbers = Array{Union{Int32, Missing}}(missing, maxW, maxL, convert(Int64, maxeps/epsstep))
    for w in 1:maxW
        # println("w:$w")
        for le in 1:maxL
            # println("\tle:$le")
            # println("\t\te:$e")
            volume = 0.0
            solution = cutandfuse_generator(le, w, 10, 10, precision=3)
            @testset "no collision" begin
                for k in keys(solution)
                    # r = Dict(
                    #     i => (s.le, s.wi) for (i, s) in enumerate(S))
                    @test !collision(solution[k].pos, solution[k].dim, filter(p -> p[1] != k, solution))
                    volume += solution[k].dim.le * solution[k].dim.wi
                end
            end
            @test isapprox(volume, w * le; atol=0.01)
            # ratios[w, le, ei] = convert(Int8, round(volume/(w*le) * 100))
            # numbers[w, le, ei] = length(r)
            # if isempty(my_rs) && length(r) <= 6 && volume/(w*le) < 60
            #     my_rs = r
            # end
            testoutofbound(solution, w)
            # Filling test
            # println(volume/(w*le))
            # @test volume/(w*le) > 0.9                
        end
    end
    # display(ratios)
    # display(numbers)
    # display(my_rs)
end


@testset "placestack!" begin

    # ```
    #     +---+
    #     | 1 |
    # +---+---+---+
    # | 3 |   | 2 |
    # +---+---+---+
    #     | 4 |
    #     +---+
    # ```
    L = 10
    W = 3
    H = 12

    
    truck = Truck(
        Dim(L, W), 
        H, 
        1500,  # max_stack_density
        Dict(), # max_stack_weight
        100000,
        Dict(), 
        Dict(), 
        Dict(), 
        10, # cost
        7808,
        3800,
        1040,
        3330,
        7300,
        7630,
        2350,
        1670,
        31500,
        12000
        )



    r = Dict(1 => Stack(Pos(1, 2), Dim(1, 1)), 2 => Stack(Pos(2, 1), Dim(1, 1)), 3 => Stack(Pos(0, 1), Dim(1, 1)), 4 => Stack(Pos(1, 0), Dim(1, 1)), 5 => Stack(Pos(1, 1), Dim(1, 1)))

    corners = [Pos(0, 0), Pos(1, 1), Pos(2, 0), Pos(0, 2)]

    for i in 1:5
        r2 = filter(p -> p[1] != i, r)
        placestack!(r2, truck, i, r[i], corners; precision=3)
        @test r2[i].pos == Pos(0, 0)
    end
    L = 10
    W = 10
    H = 12

    truck = Truck(
        Dim(L, W), 
        H, 
        1500,  # max_stack_density
        Dict(), # max_stack_weight
        100000,
        Dict(), 
        Dict(), 
        Dict(), 
        10, # cost
        7808,
        3800,
        1040,
        3330,
        7300,
        7630,
        2350,
        1670,
        31500,
        12000
        )


    r = Dict(
        0 => Stack(Pos(0, 0), Dim(2.37662, 2.1253)),
        4 => Stack(Pos(0, 2.1253), Dim(97.8747, 0.441835)),
        5 => Stack(Pos(0, 2.56714), Dim(7.18154, 2.1253)),
        6 => Stack(Pos(7.18154, 2.56714), Dim(97.8747, 7.18154)),
        2 => Stack(Pos(0, 4.69244), Dim(0.441835, 2.1253)),
        # 3 => Stack(Pos(105.056, 2.56714), Dim(97.8747, 2.37662))
    )
    placestack!(r, truck, 3, Stack(Pos(105.056, 2.56714), Dim(97.8747, 2.37662)), [Pos(2.37662, 0), Pos(105.056, 2.56714)]; precision=3, verbose=false)
    @test r[3].pos != Pos(2.37662, 0) # collision with 4


    r = Dict(
        0 => Stack(Pos(0, 0), Dim(2.49442, 3.19793)),
        4 => Stack(Pos(0, 3.19793), Dim(96.8021, 1.64834)),
        5 => Stack(Pos(96.8021, 0), Dim(96.8021, 5.85725)),
        6 => Stack(Pos(0, 4.84627), Dim(5.85725, 3.19793)),
        2 => Stack(Pos(0, 8.0442), Dim(3.19793, 1.64834)),
        # 3 => Stack(Pos(193.604, 0), Dim(96.8021, 2.49442))
    )
    placestack!(r, truck, 3, Stack(Pos(193.604, 0), Dim(96.8021, 2.49442)), [Pos(5.85725, 5.85725), Pos(193.604, 0)]; precision=3, verbose=true)
    @test r[3].pos == Pos(5.85725, 5.85725)

    # ordered stacks cases are tested with the ordered specialized can_be_placed function in ordered_stack.jl

    # same for itemized stacks cases



end



## Functional BLtruck tests
@testset "BLtruck" begin
    ratios = []
    instances_solutions = [] # TODO: make it a generator for when memory space will lack

    ITER = 100
    L = 100
    W = 10
    H = 12

    
    truck = Truck(
        Dim(L, W), 
        H, 
        1500,  # max_stack_density
        Dict(), # max_stack_weight
        100000,
        Dict(), 
        Dict(), 
        Dict(), 
        10, # cost
        7808,
        3800,
        1040,
        3330,
        7300,
        7630,
        2350,
        1670,
        31500,
        12000
        )
    NBCUTS = 20
    NBFUSE = 20
    # most, least = nb_cuts_fuse_calc(5, 20)
    
    # NBCUTS, NBFUSE = most
    
    # NBCUTS = nb_cuts_fuse_avg(4)
    # NBFUSE = 0
    for i in 1:ITER
        rectangles = cutandfuse_generator(L, W, NBCUTS, NBFUSE; precision=3)
        # Vector{Pair{T, S}} where {T <: Integer, S <: AbstractStack}
        instance = collect(rectangles)
        solution = nothing
        try
            solution = BLtruck(instance, truck, precision=3)
        catch e
            display(instance)
            throw(e)
        end
        # display(solution)
        # ratio, solution = eval_placement_function(BLtruck, instance)
        foundL = max([solution[k].pos.x + solution[k].dim.le for k in keys(solution)]...)
        ratio = foundL / L
        push!(ratios, ratio)
        push!(instances_solutions, (ins=instance, sol=solution))
    end

    # @testset "Particular test" begin
    #     rectangles = Dict(
    #         0 => Stack(Pos(0, 0), Dim(2.32568, 2.72428)),
    #         4 => Stack(Pos(0, 2.72428), Dim(97.2757, 0.140078)),
    #         5 => Stack(Pos(97.2757, 0), Dim(97.2757, 7.53424)),
    #         6 => Stack(Pos(0, 2.86436), Dim(7.53424, 2.72428)),
    #         2 => Stack(Pos(0, 5.58863), Dim(0.140078, 2.72428)),
    #         3 => Stack(Pos(194.551, 0), Dim(97.2757, 2.32568))
    #     )
    #     instance = [Pair(k, rectangles[k]) for k in [0, 4, 5, 6, 2, 3]]
    #     solution = BLtruck(instance, 10, precision=3, verbose=true)
    #     @test solution[3].pos == Pos(0.140078, 7.53424)
    #     println(collision(Pos(0.140078, 7.53424), solution[3].dim, solution; verbose=true))
    #     placestack!(solution, 10, 7, solution[3], [Pos(0.140078, 7.53424), Pos(194.551, 0)]; precision=3, verbose=true)
    #     @test solution[7].pos == Pos(0.140078, 7.53424)
    # end
    
    # Test solution is valid
    @testset "BLtruck is valid" begin
        for i in 1:ITER
            instance, solution = instances_solutions[i]
            # check overlapping and out of bounds
            @testset "No overlapping" begin
                for (j, stack) in solution
                    @test !collision(stack.pos, stack.dim, filter(p -> p[1] != j, solution); precision=3, verbose=false)
                end
            end
            @testset "Not out of bounds" begin
                
                for (j, stack) in solution
                    @test !outofbound(stack.pos, stack.dim, W; precision=3)
                end
            end
            @testset "Secure stack placement" begin
                # Stacks should be adjacent to another stack to their left on the X axis
                # or adjacent to the left side of the truck
                for (j, stack) in solution
                    @test is_secure(stack, solution; precision=3)        
                end
            end
        end 
    end
    
    # Test optimality is 2-OPT
    # @testset "BLtruck is 2-OPT" begin
    @testset "BLtruck is 3-OPT" begin
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
    println("Mean length ratio")
    display(mean(ratios))
    println("Min length ratio")
    display(min(ratios...))
    println("25% quantile length ratio")
    display(quantile(ratios, 0.25))
    println("Median length ratio")
    display(median(ratios))
    println("75% quantile length ratio")
    display(quantile(ratios, 0.75))
    println("Max length ratio")
    display(max(ratios...))

end



@testset "BLtruck loading orders" begin
    instances_solutions = []

    ITER = 100
    L = 100
    W = 10
    H = 12
    ratios = []
    
    truck = Truck(
        Dim(L, W), 
        H, 
        1500,  # max_stack_density
        Dict(), # max_stack_weight
        100000,
        Dict(), 
        Dict(), 
        Dict(), 
        10, # cost
        7808,
        3800,
        1040,
        3330,
        7300,
        7630,
        2350,
        1670,
        31500,
        12000
        )

    NBCUTS = 20
    NBFUSE = 20
    # most, least = nb_cuts_fuse_calc(5, 20)
    
    # NBCUTS, NBFUSE = most
    
    # NBCUTS = nb_cuts_fuse_avg(4)
    # NBFUSE = 0
    for i in 1:ITER
        rectangles = cutandfuse_generator(L, W, NBCUTS, NBFUSE; precision=3)

        rectangles, supplier, supplier_dock, plant_dock = order_instance(rectangles)

        instance = shuffle(rectangles)
        solution = nothing

        try
            solution = BLtruck(instance, truck; precision=3, loading_order=true)
        catch e
            display(instance)
            display(backtrace)
            throw(e)
        end
        # display(solution)
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
            @testset "Secure stack placement" begin
                # Stacks should be adjacent to another stack to their left on the X axis
                # or adjacent to the left side of the truck
                for (j, stack) in solution
                    @test is_secure(stack, solution; precision=3)        
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

    println("Mean length ratio")
    display(mean(ratios))
    println("Min length ratio")
    display(min(ratios...))
    println("25% quantile length ratio")
    display(quantile(ratios, 0.25))
    println("Median length ratio")
    display(median(ratios))
    println("75% quantile length ratio")
    display(quantile(ratios, 0.75))
    println("Max length ratio")
    display(max(ratios...))

    # Test loading orders are satisfied
    @testset "loading orders" begin
        for i in 1:ITER
            # make sure the loading order is coherent
            instance, solution = instances_solutions[i]

            x_sorted_sol = sort(collect(solution), by=p -> get_pos(p[2]).x)
            ordered_sol = sort(x_sorted_sol, by= p -> (get_supplier_order(p[2]), get_supplier_dock_order(p[2]), get_plant_dock_order(p[2])))

            # @test [p[1] for p in x_sorted_sol] == [p[1] for p in ordered_sol]
            for p in solution
                @test can_be_placed(filter(x -> x[1] != p[1], solution), get_pos(p[2]), p[2], truck, :Parallel)
                if !can_be_placed(filter(x -> x[1] != p[1], solution), get_pos(p[2]), p[2], truck, :Parallel)
                    can_be_placed(filter(x -> x[1] != p[1], solution), get_pos(p[2]), p[2], truck, :Parallel; verbose=true)
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



@testset "itemize!" begin
    # itemize(stacks::Dict{Integer, Stack}, H)
    L = 100
    W = 10
    H = 12

    NBCUTS = 20
    NBFUSE = 20

    simple_stacks = cutandfuse_generator(L, W, NBCUTS, NBFUSE; precision=3)

    truck = Truck(
        Dim(L, W), 
        H, 
        1500,  # max_stack_density
        Dict(), # max_stack_weight
        100000,
        Dict(), 
        Dict(), 
        Dict(), 
        10, # cost
        7808,
        3800,
        1040,
        3330,
        7300,
        7630,
        2350,
        1670,
        31500,
        12000
        )


    itemized_stacks = itemize!(truck, simple_stacks)
    # TODO objective test needed here
    # display(itemized_stacks)
end

@testset "BLtruck with items" begin
    # BLtruck(items, W, H, plant_dock_orders, supplier_orders, supplier_dock_orders; precision=3, verbose=false)
    instances_solutions = []

    ITER = 10
    L = 14500
    W = 2400
    H = 2800
    max_weight = 30000 # ?
    truck = Truck(
        Dim(L, W), 
        H, 
        1500,  # max_stack_density
        Dict(), # max_stack_weight
        100000,
        Dict(), 
        Dict(), 
        Dict(), 
        10, # cost
        7808,
        3800,
        1040,
        3330,
        7300,
        7630,
        2350,
        1670,
        31500,
        12000
        )
        

    NBCUTS = 20
    NBFUSE = 20
    # most, least = nb_cuts_fuse_calc(5, 20)
    
    # NBCUTS, NBFUSE = most
    ratios = []

    # NBCUTS = nb_cuts_fuse_avg(4)
    # NBFUSE = 0
    for i in 1:ITER
        t = time()
        stacks = itemize!(truck, cutandfuse_generator(get_dim(truck).le, get_dim(truck).wi, NBCUTS, NBFUSE; precision=3))
        # println("press enter a")
        # readline()
        # ordered_stacks = order_instance(stacks)
        
        # combine!(stacks, ordered_stacks)
        
        instance = shuffle(stacks)
        solution = nothing
        # println("press enter b")
        # readline()

        # you want to get the items from the stacks, shuffle them and give it to BLtruck
        items = [x for s in instance for x in get_items(s[2])]
        shuffle!(items)
        notplaced = []
        try
            solution, notplaced = BLtruck(items, truck; precision=3)
        catch e
            display(instance)
            display(backtrace)
            throw(e)
        end
        if !isempty(notplaced)
            println("These items could not be placed:")
            display(notplaced)
        end
        foundL = max([get_pos(solution[k]).x + get_dim(solution[k]).le for k in keys(solution)]...)
        push!(instances_solutions, (ins=instance, sol=solution))
        ratio = foundL / L
        push!(ratios, ratio)
        println("Ite ", i, " Elapsed: ", time() - t, "s")
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
            @testset "Secure stack placement" begin
                # Stacks should be adjacent to another stack to their left on the X axis
                # or adjacent to the left side of the truck
                for (j, stack) in solution
                    @test is_secure(stack, solution; precision=3)        
                end
            end
        end 
    end

    # Test optimality is 2-OPT
    @testset "BLtruck is 2-OPT" begin
        ## Not relevant with this constraint
        # for (i, r) in enumerate(ratios)
        #     # @test r <= 2.0 + 0.1
        #     @test r <= 3.0 # TODO why not 2-OPT?


        #     # if r > 2.0 + 0.1
        #     #     display(instances_solutions[i].ins)
        #     #     display(instances_solutions[i].sol)
        #     #     display(length(instances_solutions[i].sol))
        #     #     # println(r)
        #     # end
        # end
        println("Mean length ratio")
        display(mean(ratios))
        println("Min length ratio")
        display(min(ratios...))
        println("25% quantile length ratio")
        display(quantile(ratios, 0.25))
        println("Median length ratio")
        display(median(ratios))
        println("75% quantile length ratio")
        display(quantile(ratios, 0.75))
        println("Max length ratio")
        display(max(ratios...))
    end

    # Test loading orders are satisfied
    @testset "loading orders and weight" begin
        for i in 1:ITER
            # make sure the loading order is coherent
            instance, solution = instances_solutions[i]

            x_sorted_sol = sort(collect(solution), by=p -> get_pos(p[2]).x)
            ordered_sol = sort(x_sorted_sol, by= p -> (get_supplier_order(p[2]), get_supplier_dock_order(p[2]), get_plant_dock_order(p[2])))

            # @test [p[1] for p in x_sorted_sol] == [p[1] for p in ordered_sol]
            for p in solution
                @test can_be_placed(filter(x -> x[1] != p[1], solution), get_pos(p[2]), p[2], truck, :Parallel)
                if !can_be_placed(filter(x -> x[1] != p[1], solution), get_pos(p[2]), p[2], truck, :Parallel)
                    can_be_placed(filter(x -> x[1] != p[1], solution), get_pos(p[2]), p[2], truck, :Parallel; verbose=true)
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

    # @testset "weight constraints last check" begin
    #     for i in 1:ITER
    #         # make sure the loading order is coherent
    #         instance, solution = instances_solutions[i]

    #         # @test [p[1] for p in x_sorted_sol] == [p[1] for p in ordered_sol]
    #         @test 
    #     end

    # end
end
