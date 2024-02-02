using Test

include("../src/to_csv.jl")
include("../src/dim.jl")
include("../src/truck.jl")
include("../src/placement_algorithms.jl")


@testset "solution_to_csv" begin
    ITER = 10
    L = 14500
    W = 2400
    H = 2800
    max_weight = 30000 # ?
    truck = Truck(
        "mytruck",
        Dim(L, W), 
        H, 
        1500,  # max_stack_density
        100000, # max_stack_weight
        Dict(), 
        Dict(), 
        Dict(), 
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
    # TODO set id to stacks
    # TODO set id to items
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

    ratio = foundL / L
    println("Elapsed: ", time() - t, "s")

    directory = "test_visu"
    if !isdir(directory)
        mkdir(directory)
    end
    # reassign ids because BLtruck creates new stacks id-less
    for (i, stack) in solution
        set_id!(stack, string(get_id(truck), "_", i))
    end

    solution_to_csv(truck, solution, directory, append=false)
end