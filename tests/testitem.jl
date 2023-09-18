using Test
using Random

include("../src/item.jl")

@testset "rand_items" begin
    n = 100
    min_products = 1
    max_products = 10
    max_h = 100
    max_w = 100
    max_items_per_stack = n
    L = 100
    W = 10
    plant = randstring(8)

    products = rand_products(min_products, max_products, max_w * max_items_per_stack, max_items_per_stack)

    items = rand_items(
        n, 
        products,
        max_h, 
        max_w, 
        L, 
        W, 
        plant; 
        min_dim=0.001
        )

    # display(items)

end

@testset "add_item!" begin
    n = 100
    min_products = 1
    max_products = 10
    max_h = 100
    max_w = 100
    max_items_per_stack = n
    L = 100
    W = 10
    plant = randstring(8)
    products = rand_products(min_products, max_products, max_w * max_items_per_stack, max_items_per_stack)

    items = rand_items(
        n, 
        products,
        max_h, 
        max_w, 
        L, 
        W, 
        plant; 
        min_dim=0.001
        )

    newstack = ItemizedStack(
        1, 
        2,
        3
        )

    display(newstack)
    add_item!(newstack, items[1])
    display(newstack)
    add_item!(newstack, items[2])
    display(newstack)

    # readline()
end

@testset "valid_stack" begin

    time_window = (earliest=0, latest=10)
    dim=Dim(1, 1)
    # pos::Pos
    height=1.0
    weight=1.0
    stackability_code="A"
    forced_orientation=[:Free, :Horizontal, :Vertical]
    plant="PA"
    plant_dock="PDA"
    supplier="SA"
    supplier_dock="SDA"
    inventory_cost=1.0
    nesting_height=0.0
    product = Product(10, 3.0)

    i = Item(
        time_window,
        dim,
        # pos::Pos,
        height,
        weight,
        stackability_code,
        forced_orientation[1],
        plant,
        plant_dock,
        supplier,
        supplier_dock,
        inventory_cost,
        nesting_height,
        product
        )
    stack = ItemizedStack(1, 2, 3)
    add_item!(stack, i)

    # test height
    @test valid_stack(stack, i, 10)
    @test !valid_stack(stack, i, 1)

    add_item!(stack, copy(i))
    # test weight
    
    
    i2 = Item(
        time_window,
        dim,
        # pos::Pos,
        height,
        weight*2,
        stackability_code,
        forced_orientation[1],
        plant,
        plant_dock,
        supplier,
        supplier_dock,
        inventory_cost,
        nesting_height,
        product
        )
    @test !valid_stack(stack, i2, 10)

    add_item!(stack, copy(i))

    # test number
    @test !valid_stack(stack, copy(i), 10)
    display(stack)

    # TODO test orientation
    stack = ItemizedStack(1, 2, 3)
    add_item!(stack, i)

    i = Item(
        time_window,
        dim,
        # pos::Pos,
        height,
        weight,
        stackability_code,
        forced_orientation[2],
        plant,
        plant_dock,
        supplier,
        supplier_dock,
        inventory_cost,
        nesting_height,
        product
        )
    # stack free, item horizontal
    @test valid_stack(stack, copy(i), 10)


    stack = ItemizedStack(1, 2, 3)
    add_item!(stack, i)

    # stack horizontal, item horizontal
    @test valid_stack(stack, copy(i), 10)

    i = Item(
        time_window,
        dim,
        # pos::Pos,
        height,
        weight,
        stackability_code,
        forced_orientation[1],
        plant,
        plant_dock,
        supplier,
        supplier_dock,
        inventory_cost,
        nesting_height,
        product
        )
    
    # stack horizontal, item free
    @test valid_stack(stack, copy(i), 10)

    i = Item(
        time_window,
        dim,
        # pos::Pos,
        height,
        weight,
        stackability_code,
        forced_orientation[3],
        plant,
        plant_dock,
        supplier,
        supplier_dock,
        inventory_cost,
        nesting_height,
        product
        )
        # stack horizontal, item vertical
    @test !valid_stack(stack, copy(i), 10)

end

@testset "is_candidate_stack" begin
    time_window = (earliest=0, latest=10)
    dim=Dim(1, 1)
    # pos::Pos
    height=1.0
    weight=1.0
    stackability_code=["A", "B"]
    forced_orientation=:Free
    plant="PA"
    plant_dock=["PDA", "PDB"]
    supplier=["SA", "SB"]
    supplier_dock=["SDA", "SDB"]
    inventory_cost=1.0
    nesting_height=0.0
    product = Product(10, 3.0)

    i1 = Item(
        time_window,
        dim,
        # pos::Pos,
        height,
        weight,
        stackability_code[1],
        forced_orientation,
        plant,
        plant_dock[1],
        supplier[1],
        supplier_dock[1],
        inventory_cost,
        nesting_height,
        product
        )
    
    items = [Item(
        time_window,
        dim,
        # pos::Pos,
        height,
        weight,
        stackability_code[a],
        forced_orientation,
        plant,
        plant_dock[b],
        supplier[c],
        supplier_dock[d],
        inventory_cost,
        nesting_height,
        product
        ) for (a, b, c, d) in [(1, 1, 1, 1), (1, 1, 1, 2), (1, 1, 2, 1), (1, 2, 1, 1), (2, 1, 1, 1)]]
    

    # all good
    stack = ItemizedStack(1, 2, 3)
    add_item!(stack, items[1])
    @test is_candidate_stack(stack, i1)

    # supplier dock
    stack = ItemizedStack(1, 2, 3)
    add_item!(stack, items[2])
    @test !is_candidate_stack(stack, i1)

    # supplier
    stack = ItemizedStack(1, 2, 3)
    add_item!(stack, items[3])
    @test !is_candidate_stack(stack, i1)

    # plant dock
    stack = ItemizedStack(1, 2, 3)
    add_item!(stack, items[4])
    @test !is_candidate_stack(stack, i1)

    # stackability code
    stack = ItemizedStack(1, 2, 3)
    add_item!(stack, items[5])
    @test !is_candidate_stack(stack, i1)


end

@testset "adding weight to stack" begin
    mystack = ItemizedStack(1, 1, 1)

    item = Item(
    (earliest=1, latest=10),
    Dim(1, 1),
    # pos::Pos,
    1,
    10, # weight
    "",
    :Free,
    "",
    "",
    "",
    "",
    1,
    0,
    Product(1000, 1000)
    )

    @test get_weight(mystack) == 0
    
    add_item!(mystack, copy(item))
    
    @test get_weight(mystack) == 10
    
    
    add_item!(mystack, copy(item))

    @test get_weight(mystack) == 20
    
    add_item!(mystack, copy(item))
    
    @test get_weight(mystack) == 30
    
    add_item!(mystack, copy(item))
    
    @test get_weight(mystack) == 40
end

@testset "make_stacks1" begin

    n = 100
    min_products = 1
    max_products = 10
    max_h = 100
    max_w = 100
    max_items_per_stack = n
    L = 100
    W = 10
    plant = randstring(8)
    products = rand_products(min_products, max_products, max_w * max_items_per_stack, max_items_per_stack)

    items = rand_items(
        n, 
        products, 
        max_h, 
        max_w, 
        L, 
        W, 
        plant; 
        min_dim=0.001
        )

    # find plant docks and assign random order
    plant_dock_orders = Dict()
    supplier_orders = Dict() # key is supplier
    supplier_dock_orders = Dict() # key is supplier then supplier_dock
    for i in items
        if !(get_plant_dock(i) in keys(plant_dock_orders))
            plant_dock_orders[get_plant_dock(i)] = nothing
        end
        if !(get_supplier(i) in keys(supplier_orders))
            supplier_orders[get_supplier(i)] = nothing
            supplier_dock_orders[get_supplier(i)] = Dict()
        end
        if !(get_supplier_dock(i) in keys(supplier_dock_orders[get_supplier(i)]))
            supplier_dock_orders[get_supplier(i)][get_supplier_dock(i)] = nothing
        end
    end

    orders = collect(1:length(keys(plant_dock_orders)))
    shuffle!(orders)
    for plant_dock in keys(plant_dock_orders)
        plant_dock_orders[plant_dock] = pop!(orders)
    end

    # find supplier and assign random order
    orders = collect(1:length(keys(supplier_orders)))
    shuffle!(orders)
    # for each supplier
    for supplier in keys(supplier_orders)
        supplier_orders[supplier] = pop!(orders)
       
        dock_orders = collect(1:length(keys(supplier_dock_orders[supplier])))
        shuffle!(dock_orders)

        # find supplier docks and assign random order
        for supplier_dock in keys(supplier_dock_orders[supplier])
            supplier_dock_orders[supplier][supplier_dock] = pop!(dock_orders)
        end
    end

    # display(plant_dock_orders)
    # display(supplier_orders)
    # display(supplier_dock_orders)
    
    max_height = 100

    stacks = make_stacks(
        convert(Vector{Item}, items), 
        plant_dock_orders, 
        supplier_orders, 
        supplier_dock_orders, 
        max_height
        )
    
    for supplier in keys(stacks)
        display(stacks[supplier])
    end
    # display(stacks[collect(keys(stacks))[begin]])

    # sum items and make sure == n
    m = sum([length(get_items(s)) for supplier in keys(stacks) for s in stacks[supplier]])
    @test m == n

end


@testset "make_stacks2" begin

    n = 100
    min_products = 1
    max_products = 1
    max_h = 1
    max_w = 1
    max_items_per_stack = n
    L = 1000
    W = 1000
    plant = randstring(8)
    products = rand_products(min_products, max_products, max_w * max_items_per_stack, max_items_per_stack)

    items = rand_items(
        n, 
        products, 
        max_h, 
        max_w, 
        L, 
        W, 
        plant; 
        min_dim=0.001
        )

    # find plant docks and assign random order
    plant_dock_orders = Dict()
    supplier_orders = Dict() # key is supplier
    supplier_dock_orders = Dict() # key is supplier then supplier_dock
    for i in items
        if !(get_plant_dock(i) in keys(plant_dock_orders))
            plant_dock_orders[get_plant_dock(i)] = nothing
        end
        if !(get_supplier(i) in keys(supplier_orders))
            supplier_orders[get_supplier(i)] = nothing
            supplier_dock_orders[get_supplier(i)] = Dict()
        end
        if !(get_supplier_dock(i) in keys(supplier_dock_orders[get_supplier(i)]))
            supplier_dock_orders[get_supplier(i)][get_supplier_dock(i)] = nothing
        end
    end

    orders = collect(1:length(keys(plant_dock_orders)))
    shuffle!(orders)
    for plant_dock in keys(plant_dock_orders)
        plant_dock_orders[plant_dock] = pop!(orders)
    end

    # find supplier and assign random order
    orders = collect(1:length(keys(supplier_orders)))
    shuffle!(orders)
    # for each supplier
    for supplier in keys(supplier_orders)
        supplier_orders[supplier] = pop!(orders)
       
        dock_orders = collect(1:length(keys(supplier_dock_orders[supplier])))
        shuffle!(dock_orders)

        # find supplier docks and assign random order
        for supplier_dock in keys(supplier_dock_orders[supplier])
            supplier_dock_orders[supplier][supplier_dock] = pop!(dock_orders)
        end
    end

    # display(plant_dock_orders)
    # display(supplier_orders)
    # display(supplier_dock_orders)
    
    max_height = 100000

    stacks = make_stacks(
        convert(Vector{Item}, items), 
        plant_dock_orders, 
        supplier_orders, 
        supplier_dock_orders, 
        max_height
        )
    
    for supplier in keys(stacks)
        display(stacks[supplier])
    end
    # display(stacks[collect(keys(stacks))[begin]])

    # sum items and make sure == n
    m = sum([length(get_items(s)) for supplier in keys(stacks) for s in stacks[supplier]])
    @test m == n

end

