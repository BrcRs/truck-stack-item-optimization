using Test
using Random

include("../src/item.jl")

println("testitem.jl outputs ---")

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


@testset "dist_stacks_to_trailer with an without given stack" begin
    """
    ```
        +---+
        | 1 |
    +---+---+---+
    | 2 | 3 | 4 |
    +---+---+---+
        | 5 |
        +---+
    ```
    """

    EJ_eh = 1670
    EJ_hr = 7630

    EJ_cr = 2350
    EM = 7300
    EM_mm = 12000
    EM_mr = 31500

    CJ_fm = 3800
    CJ_fc = 1040
    CJ_fh = 3330
    CM = 7808

    _tm_t = 5

    _ej_e = (1 * (1 + (1)/2) + 1 * (0 + 1/2) + 1 * (1 + 1/2) + 1 * (2 + 1/2) + 1 * (1 + 1/2)) / _tm_t
    _ej_r = EJ_eh + EJ_hr - _ej_e

    _em_h = (_tm_t * _ej_r + EM * EJ_cr) / EJ_hr

    _em_r = _tm_t + EM - _em_h
    _em_m = (CM * CJ_fc + _em_h * CJ_fh) / CJ_fm

    truck = Truck(Dim(14500, 2400), 2800, 1500, 100000, Dict(), Dict(), Dict(), CM, CJ_fm, CJ_fc, CJ_fh, EM, EJ_hr, EJ_cr, EJ_eh, EM_mr, EM_mm)

    positions = [Pos(1, 2), Pos(0, 1), Pos(1, 1), Pos(2, 1), Pos(1, 0)]
    stacks = [ItemizedStack(1, 1, i) for i in 1:5]
    stacks = [newStack(stack, positions[i], Dim(1, 1)) for (i, stack) in enumerate(stacks)]

    max_stackability, max_weight = 10, 10
    p = Product(max_stackability, max_weight)

    items = [simpleItem(p) for i in 1:5]

    for i in 1:5
        add_item!(stacks[i], items[i])
    end

    @test dist_stacks_to_trailer(stacks, truck) == (_tm_t, _ej_e, _ej_r, _em_h, _em_r, _em_m)

    @test dist_stacks_to_trailer(stacks[1:end-1], stacks[end], truck) == (_tm_t, _ej_e, _ej_r, _em_h, _em_r, _em_m)



end

@testset "valid_axle_pressure" begin
    # valid_axle_pressure(stacks, s::ItemizedStack, truck::Truck; fastexit=false, precision=3)
    """
    ```
        +---+
        | 1 |
    +---+---+---+
    | 2 | 3 | 4 |
    +---+---+---+
        | 5 |
        +---+
    ```
    """

    EJ_eh = 1670
    EJ_hr = 7630

    EJ_cr = 2350
    EM = 7300
    EM_mm = 12000
    EM_mr = 31500

    CJ_fm = 3800
    CJ_fc = 1040
    CJ_fh = 3330
    CM = 7808

    _tm_t = 5

    # tmp = ((1 + (1)/2) + (0 + 1/2) + (1 + 1/2) + (2 + 1/2) + (1 + 1/2)) * (_tm_t / 5)/ _tm_t
    # tmp = ((1 + (1)/2) + (0 + 1/2) + (1 + 1/2) + (2 + 1/2) + (1 + 1/2))  / 5

    # _ej_e = (1 * (1 + (1)/2) + 1 * (0 + 1/2) + 1 * (1 + 1/2) + 1 * (2 + 1/2) + 1 * (1 + 1/2)) / _tm_t
    _ej_e = ((1 + (1)/2) + (0 + 1/2) + (1 + 1/2) + (2 + 1/2) + (1 + 1/2)) / 5
    _ej_r = EJ_eh + EJ_hr - _ej_e

    _em_h = (_tm_t * _ej_r + EM * EJ_cr) / EJ_hr

    _em_r = _tm_t + EM - _em_h
    _em_m = (CM * CJ_fc + _em_h * CJ_fh) / CJ_fm

    truck = Truck(Dim(14500, 2400), 2800, 1500, 100000, Dict("A" => 1, "B" => 2), Dict(), Dict(), CM, CJ_fm, CJ_fc, CJ_fh, EM, EJ_hr, EJ_cr, EJ_eh, EM_mr, EM_mm)

    positions = [Pos(1, 2), Pos(0, 1), Pos(1, 1), Pos(2, 1), Pos(1, 0)]
    stacks = [ItemizedStack(i < 3 ? 1 : 2, 1, 1) for i in 1:5]
    stacks = [newStack(stack, positions[i], Dim(1, 1)) for (i, stack) in enumerate(stacks)]

    max_stackability, max_weight = 10, 10
    p = Product(max_stackability, max_weight)

    items = [simpleItem(p, supplier=i < 3 ? "A" : "B") for i in 1:5]

    for i in 1:5
        add_item!(stacks[i], items[i])
    end

    # test true and false

    # Change EM_mm and EM_mr to test false and true cases
    
    # true
    @test valid_axle_pressure(stacks[1:end-1], stacks[end], truck; fastexit=false, precision=3)

    # false for supplier A alone
    # we want em_m > EM_mm
    # (CM * CJ_fc + _em_h * CJ_fh) / CJ_fm > EM_mm
    # _em_h * CJ_fh / CJ_fm > EM_mm + CM * CJ_fc / CJ_fm
    # _em_h > (CJ_fm * EM_mm + CM * CJ_fc) / CJ_fh
    # (_tm_tA * _ej_r + EM * EJ_cr) / EJ_hr > (CJ_fm * EM_mm + CM * CJ_fc) / CJ_fh
    # _tm_tA * _ej_r + EM * EJ_cr > (CJ_fm * EM_mm + CM * CJ_fc) * EJ_hr / CJ_fh
    # _tm_tA * _ej_r > (CJ_fm * EM_mm + CM * CJ_fc) * EJ_hr / CJ_fh - EM * EJ_cr
    # _tm_tA * _ej_r > (3800 * 12000 + 7808 * 1040) * 7630 / 3330 - 7300 * 2350
    # _tm_tA * _ej_r > 105 933 901.38138138
    # _tm_tA  > 105 933 901.38138138 / _ej_r
    # _tm_tA  > 105 933 901.38138138 / (EJ_eh + EJ_hr - _ej_e)
    # _tm_tA  > 105 933 901.38138138 / (EJ_eh + EJ_hr - (wA * (1 + (1)/2) + wA * (0 + 1/2))  / _tm_tA)
    # 2 * wA  > 105 933 901.38138138 / (EJ_eh + EJ_hr - (wA * (1 + (1)/2) + wA * (0 + 1/2))  / (2 * wA))
    # 2 * wA * (EJ_eh + EJ_hr - 2 * wA  / (2 * wA))  > 105 933 901.38138138
    # 2 * wA * (EJ_eh + EJ_hr - 1)  > 105 933 901.38138138
    # 2 * wA  > 105 933 901.38138138 / (EJ_eh + EJ_hr - 1)
    # 2 * wA  > 105 933 901.38138138 / (1670 + 7629)
    # wA  > 105 933 901.38138138 / (1670 + 7629) * 2
    # wA  > 22783.934053421093

    
    EJ_eh = 1670
    EJ_hr = 7630

    EJ_cr = 2350
    EM = 7300
    EM_mm = 12000
    EM_mr = 31500

    CJ_fm = 3800
    CJ_fc = 1040
    CJ_fh = 3330
    CM = 7808

    # _tm_t = 13000

    # tmp = ((1 + (1)/2) + (0 + 1/2) + (1 + 1/2) + (2 + 1/2) + (1 + 1/2)) * (_tm_t / 5)/ _tm_t
    # tmp = ((1 + (1)/2) + (0 + 1/2) + (1 + 1/2) + (2 + 1/2) + (1 + 1/2))  / 5

    # _ej_e = (1 * (1 + (1)/2) + 1 * (0 + 1/2) + 1 * (1 + 1/2) + 1 * (2 + 1/2) + 1 * (1 + 1/2)) / _tm_t
    # _ej_e = ((1 + (1)/2) + (0 + 1/2) + (1 + 1/2) + (2 + 1/2) + (1 + 1/2)) / 5
    # _ej_r = EJ_eh + EJ_hr - _ej_e

    # _em_h = (_tm_t * _ej_r + EM * EJ_cr) / EJ_hr

    # _em_r = _tm_t + EM - _em_h
    # _em_m = (CM * CJ_fc + _em_h * CJ_fh) / CJ_fm

    truck = Truck(Dim(14500, 2400), 2800, 1500, 100000, Dict("A" => 1, "B" => 2), Dict(), Dict(), CM, CJ_fm, CJ_fc, CJ_fh, EM, EJ_hr, EJ_cr, EJ_eh, EM_mr, EM_mm)

    positions = [Pos(1, 2), Pos(0, 1), Pos(1, 1), Pos(2, 1), Pos(1, 0)]
    stacks = [ItemizedStack(i < 3 ? 1 : 2, 1, 1) for i in 1:5]
    stacks = [newStack(stack, positions[i], Dim(1, 1)) for (i, stack) in enumerate(stacks)]

    max_stackability, max_weight = 10, 10
    p = Product(max_stackability, max_weight)

    items = [simpleItem(p, supplier=i < 3 ? "A" : "B", weight=supplier=i < 3 ? 22784 : 1) for i in 1:5]

    for i in 1:5
        add_item!(stacks[i], items[i])
    end

    @test !valid_axle_pressure(stacks[1:end-1], stacks[end], truck; fastexit=false, precision=3, verbose=true)
    
    
    
    # false for both
    # we want em_m > EM_mm
    # (CM * CJ_fc + _em_h * CJ_fh) / CJ_fm > EM_mm
    # _em_h * CJ_fh / CJ_fm > EM_mm + CM * CJ_fc / CJ_fm
    # _em_h > (CJ_fm * EM_mm + CM * CJ_fc) / CJ_fh
    # (_tm_t * _ej_r + EM * EJ_cr) / EJ_hr > (CJ_fm * EM_mm + CM * CJ_fc) / CJ_fh
    # _tm_t * _ej_r + EM * EJ_cr > (CJ_fm * EM_mm + CM * CJ_fc) * EJ_hr / CJ_fh
    # _tm_t * _ej_r > (CJ_fm * EM_mm + CM * CJ_fc) * EJ_hr / CJ_fh - EM * EJ_cr
    # _tm_t * _ej_r > (3800 * 12000 + 7808 * 1040) * 7630 / 3330 - 7300 * 2350
    # _tm_t * _ej_r > 105 933 901.38138138
    # _tm_t  > 105 933 901.38138138 / _ej_r
    # _tm_t  > 105 933 901.38138138 / (EJ_eh + EJ_hr - _ej_e)
    # _tm_t  > 105 933 901.38138138 / (EJ_eh + EJ_hr - ((1 + (1)/2) + (0 + 1/2) + (1 + 1/2) + (2 + 1/2) + (1 + 1/2))  / 5)
    # _tm_t  > 105 933 901.38138138 / (1670 + 7630 - ((1 + (1)/2) + (0 + 1/2) + (1 + 1/2) + (2 + 1/2) + (1 + 1/2))  / 5)
    # _tm_t  > 11392.579596857706
    # weight per item = 2278.5159


    EJ_eh = 1670
    EJ_hr = 7630

    EJ_cr = 2350
    EM = 7300
    EM_mm = 12000
    EM_mr = 31500

    CJ_fm = 3800
    CJ_fc = 1040
    CJ_fh = 3330
    CM = 7808

    _tm_t = 13000

    # tmp = ((1 + (1)/2) + (0 + 1/2) + (1 + 1/2) + (2 + 1/2) + (1 + 1/2)) * (_tm_t / 5)/ _tm_t
    # tmp = ((1 + (1)/2) + (0 + 1/2) + (1 + 1/2) + (2 + 1/2) + (1 + 1/2))  / 5

    # _ej_e = (1 * (1 + (1)/2) + 1 * (0 + 1/2) + 1 * (1 + 1/2) + 1 * (2 + 1/2) + 1 * (1 + 1/2)) / _tm_t
    _ej_e = ((1 + (1)/2) + (0 + 1/2) + (1 + 1/2) + (2 + 1/2) + (1 + 1/2)) / 5
    _ej_r = EJ_eh + EJ_hr - _ej_e

    _em_h = (_tm_t * _ej_r + EM * EJ_cr) / EJ_hr

    _em_r = _tm_t + EM - _em_h
    _em_m = (CM * CJ_fc + _em_h * CJ_fh) / CJ_fm

    truck = Truck(Dim(14500, 2400), 2800, 1500, 100000, Dict("A" => 1, "B" => 2), Dict(), Dict(), CM, CJ_fm, CJ_fc, CJ_fh, EM, EJ_hr, EJ_cr, EJ_eh, EM_mr, EM_mm)

    positions = [Pos(1, 2), Pos(0, 1), Pos(1, 1), Pos(2, 1), Pos(1, 0)]
    stacks = [ItemizedStack(i < 3 ? 1 : 2, 1, 1) for i in 1:5]
    stacks = [newStack(stack, positions[i], Dim(1, 1)) for (i, stack) in enumerate(stacks)]

    max_stackability, max_weight = 10, 10
    p = Product(max_stackability, max_weight)

    items = [simpleItem(p, supplier=i < 3 ? "A" : "B", weight=_tm_t/5) for i in 1:5]

    for i in 1:5
        add_item!(stacks[i], items[i])
    end

    @test !valid_axle_pressure(stacks[1:end-1], stacks[end], truck; fastexit=false, precision=3, verbose=true)
    # truck has no suppliers, but those suppliers are looked for by valid_axle_pressure TODO


end

@testset "valid_stack" begin

    EJ_eh = 1670
    EJ_hr = 7630

    EJ_cr = 2350
    EM = 7300
    EM_mm = 12000
    EM_mr = 31500

    CJ_fm = 3800
    CJ_fc = 1040
    CJ_fh = 3330
    CM = 7808
    truck = Truck(Dim(14500, 2400), 1, 1500, 100000, Dict("A" => 1, "B" => 2), Dict(), Dict(), CM, CJ_fm, CJ_fc, CJ_fh, EM, EJ_hr, EJ_cr, EJ_eh, EM_mr, EM_mm)


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
    stack = ItemizedStack(OrderedStack(Pos(0, 0), Dim(1, 1), 1, 2, 3))
    add_item!(stack, i)

    # test height
    truck = set_height(truck, 10)
    @test valid_stack([], stack, i, truck; verbose=true)
    truck = set_height(truck, 1)
    @test !valid_stack([], stack, i,truck)

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

    truck = set_height(truck, 10)
    @test !valid_stack([], stack, i2, truck)

    add_item!(stack, copy(i))

    # test number
    @test !valid_stack([], stack, copy(i), truck)
    display(stack)

    # TODO test orientation
    stack = ItemizedStack(OrderedStack(Pos(0, 0), Dim(1, 1), 1, 2, 3))
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
    @test valid_stack([], stack, copy(i), truck)


    stack = ItemizedStack(OrderedStack(Pos(0, 0), Dim(1, 1), 1, 2, 3))
    add_item!(stack, i)

    # stack horizontal, item horizontal
    @test valid_stack([], stack, copy(i), truck)

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
    @test valid_stack([], stack, copy(i), truck)

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
    @test !valid_stack([], stack, copy(i), truck)

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

# @testset "make_stacks1" begin

#     n = 100
#     min_products = 1
#     max_products = 10
#     max_h = 100
#     max_w = 100
#     max_items_per_stack = n
#     L = 100
#     W = 10
#     plant = randstring(8)
#     products = rand_products(min_products, max_products, max_w * max_items_per_stack, max_items_per_stack)

#     items = rand_items(
#         n, 
#         products, 
#         max_h, 
#         max_w, 
#         L, 
#         W, 
#         plant; 
#         min_dim=0.001
#         )

#     # find plant docks and assign random order
#     plant_dock_orders = Dict()
#     supplier_orders = Dict() # key is supplier
#     supplier_dock_orders = Dict() # key is supplier then supplier_dock
#     for i in items
#         if !(get_plant_dock(i) in keys(plant_dock_orders))
#             plant_dock_orders[get_plant_dock(i)] = nothing
#         end
#         if !(get_supplier(i) in keys(supplier_orders))
#             supplier_orders[get_supplier(i)] = nothing
#             supplier_dock_orders[get_supplier(i)] = Dict()
#         end
#         if !(get_supplier_dock(i) in keys(supplier_dock_orders[get_supplier(i)]))
#             supplier_dock_orders[get_supplier(i)][get_supplier_dock(i)] = nothing
#         end
#     end

#     orders = collect(1:length(keys(plant_dock_orders)))
#     shuffle!(orders)
#     for plant_dock in keys(plant_dock_orders)
#         plant_dock_orders[plant_dock] = pop!(orders)
#     end

#     # find supplier and assign random order
#     orders = collect(1:length(keys(supplier_orders)))
#     shuffle!(orders)
#     # for each supplier
#     for supplier in keys(supplier_orders)
#         supplier_orders[supplier] = pop!(orders)
       
#         dock_orders = collect(1:length(keys(supplier_dock_orders[supplier])))
#         shuffle!(dock_orders)

#         # find supplier docks and assign random order
#         for supplier_dock in keys(supplier_dock_orders[supplier])
#             supplier_dock_orders[supplier][supplier_dock] = pop!(dock_orders)
#         end
#     end

#     # display(plant_dock_orders)
#     # display(supplier_orders)
#     # display(supplier_dock_orders)
    
#     max_height = 100

#     stacks = make_stacks(
#         convert(Vector{Item}, items), 
#         plant_dock_orders, 
#         supplier_orders, 
#         supplier_dock_orders, 
#         max_height
#         )
    
#     for supplier in keys(stacks)
#         display(stacks[supplier])
#     end
#     # display(stacks[collect(keys(stacks))[begin]])

#     # sum items and make sure == n
#     m = sum([length(get_items(s)) for supplier in keys(stacks) for s in stacks[supplier]])
#     @test m == n

# end


# @testset "make_stacks2" begin

#     n = 100
#     min_products = 1
#     max_products = 1
#     max_h = 1
#     max_w = 1
#     max_items_per_stack = n
#     L = 1000
#     W = 1000
#     plant = randstring(8)
#     products = rand_products(min_products, max_products, max_w * max_items_per_stack, max_items_per_stack)

#     items = rand_items(
#         n, 
#         products, 
#         max_h, 
#         max_w, 
#         L, 
#         W, 
#         plant; 
#         min_dim=0.001
#         )

#     # find plant docks and assign random order
#     plant_dock_orders = Dict()
#     supplier_orders = Dict() # key is supplier
#     supplier_dock_orders = Dict() # key is supplier then supplier_dock
#     for i in items
#         if !(get_plant_dock(i) in keys(plant_dock_orders))
#             plant_dock_orders[get_plant_dock(i)] = nothing
#         end
#         if !(get_supplier(i) in keys(supplier_orders))
#             supplier_orders[get_supplier(i)] = nothing
#             supplier_dock_orders[get_supplier(i)] = Dict()
#         end
#         if !(get_supplier_dock(i) in keys(supplier_dock_orders[get_supplier(i)]))
#             supplier_dock_orders[get_supplier(i)][get_supplier_dock(i)] = nothing
#         end
#     end

#     orders = collect(1:length(keys(plant_dock_orders)))
#     shuffle!(orders)
#     for plant_dock in keys(plant_dock_orders)
#         plant_dock_orders[plant_dock] = pop!(orders)
#     end

#     # find supplier and assign random order
#     orders = collect(1:length(keys(supplier_orders)))
#     shuffle!(orders)
#     # for each supplier
#     for supplier in keys(supplier_orders)
#         supplier_orders[supplier] = pop!(orders)
       
#         dock_orders = collect(1:length(keys(supplier_dock_orders[supplier])))
#         shuffle!(dock_orders)

#         # find supplier docks and assign random order
#         for supplier_dock in keys(supplier_dock_orders[supplier])
#             supplier_dock_orders[supplier][supplier_dock] = pop!(dock_orders)
#         end
#     end

#     # display(plant_dock_orders)
#     # display(supplier_orders)
#     # display(supplier_dock_orders)
    
#     max_height = 100000

#     stacks = make_stacks(
#         convert(Vector{Item}, items), 
#         plant_dock_orders, 
#         supplier_orders, 
#         supplier_dock_orders, 
#         max_height
#         )
    
#     for supplier in keys(stacks)
#         display(stacks[supplier])
#     end
#     # display(stacks[collect(keys(stacks))[begin]])

#     # sum items and make sure == n
#     m = sum([length(get_items(s)) for supplier in keys(stacks) for s in stacks[supplier]])
#     @test m == n

# end



# @testset "can_be_placed" begin
#     """
#     can_be_placed(solution, o::Pos, s::ItemizedStack, truck::Truck, orientation::Symbol; precision=3, verbose=false)

#     Return true if stack `s` can be placed in partial solution `solution` without breaking any constraint.
#     """
#     # can_be_placed(solution, o::Pos, s::ItemizedStack, truck::Truck, orientation::Symbol; precision=3, verbose=false)

#     # make simple stacks to add to solution
#     stacks = Dict(
#         i => ItemizedStack(1, 1, i) for i in 1:5
#         )
#     max_stackability, max_weight = 10, 10
#     p = Product(1000, 1000)

#     items = [simpleItem()]

#     add_item!(stack, i)
#     # find corners
#     error("WIP")
#     # make a test stack

#     # make a truck

#     # determine orientation

#     # only test weight related constraints

#     can_be_placed(solution, o::Pos, s::ItemizedStack, truck::Truck, orientation::Symbol; precision=3, verbose=false)


# end