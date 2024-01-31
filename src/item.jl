using AutoHashEquals
using Random

include("placement.jl")
include("ordered_stacks.jl")
include("truck.jl")




@auto_hash_equals struct Product
    max_stackability::Integer # max number of items an item of product can support
    max_weight::Real # max weight an item of product can support 
end

get_max_stackability(p::Product) = p.max_stackability
get_max_weight(p::Product) = p.max_weight

@auto_hash_equals struct Item
    time_window::@NamedTuple{earliest::Any, latest::Any}
    dim::Dim
    # pos::Pos
    height::Real
    weight::Real
    stackability_code
    forced_orientation::Symbol
    plant
    plant_dock
    supplier
    supplier_dock
    inventory_cost
    nesting_height
    product::Product
end

get_time_window(i::Item) = i.time_window
get_dim(i::Item) = i.dim
get_height(i::Item) = i.height
get_weight(i::Item) = i.weight
get_stackability_code(i::Item) = i.stackability_code
get_forced_orientation(i::Item) = i.forced_orientation
get_plant(i::Item) = i.plant
get_plant_dock(i::Item) = i.plant_dock
get_supplier(i::Item) = i.supplier
get_supplier_dock(i::Item) = i.supplier_dock
get_inventory_cost(i::Item) = i.inventory_cost
get_nesting_height(i::Item) = i.nesting_height
get_product(i::Item) = i.product
get_max_weight(i::Item) = get_max_weight(i.product)
get_max_stackability(i::Item) = get_max_stackability(i.product)

Base.show(io::IO, it::Item) = 
    print(io, 
        "Item(timing=", (get_time_window(it).earliest, get_time_window(it).latest), 
        ", Dim(", round(get_dim(it).le, digits=3), ", ", round(get_dim(it).wi, digits=3), 
        "), h=", round(get_height(it), digits=3), 
        ", w=", round(get_weight(it), digits=3), 
        ", stackability=", get_stackability_code(it),
        ", forient.=", get_forced_orientation(it),
        ", plant&dock=", get_plant(it), " ", get_plant_dock(it),
        # ", plant dock=", get_plant_dock(it),
        ", supp.&dock=", get_supplier(it), " ", get_supplier_dock(it),
        # ", supp. dock=", get_supplier_dock(it),
        ", cost=", get_inventory_cost(it),
        ", nesting=", get_nesting_height(it),
        ", ", get_product(it))

function Base.copy(it::Item)
    return Item(
        it.time_window,
        it.dim,
        # pos::Pos,
        it.height,
        it.weight,
        it.stackability_code,
        it.forced_orientation,
        it.plant,
        it.plant_dock,
        it.supplier,
        it.supplier_dock,
        it.inventory_cost,
        it.nesting_height,
        it.product # no deep copy necessary because immutable
        )
end

function simpleItem(product;
    time_window=(earliest=1, latest=2),
    dim=Dim(1, 1),
    # pos::Pos,
    height=1,
    weight=1,
    stackability_code="",
    forced_orientation=:Free,
    plant="",
    plant_dock="",
    supplier="",
    supplier_dock="",
    inventory_cost=1,
    nesting_height=0    
    )
    return Item(time_window,
    dim,
    height,
    weight,
    stackability_code,
    forced_orientation,
    plant,
    plant_dock,
    supplier,
    supplier_dock,
    inventory_cost,
    nesting_height,
    product)

end

"""
An ordered stack which can hold items.
"""
mutable struct ItemizedStack <: AbstractOrderedStack
    ordered_stack::Union{Nothing, OrderedStack}
    items::Vector{Item}
    weight::Real
    height::Real
    forced_orientation::Symbol
    minmax_stackability::Number
    function ItemizedStack(ordered_stack::OrderedStack,
        items::Vector{Item},
        weight::Real,
        height::Real,
        forced_orientation::Symbol,
        minmax_stackability::Number)
        if !(forced_orientation in [:Free, :Horizontal, :Vertical])
            throw(ArgumentError("forced_orientation should be taken among [:Free, :Horizontal, :Vertical].\nGot $forced_orientation."))
        end
        new(ordered_stack, items, weight, height, forced_orientation, minmax_stackability)
    end

end
 

function Base.copy(is::ItemizedStack)
    return ItemizedStack(
        is.ordered_stack,
        copy(is.items),
        is.weight,
        is.height,
        is.forced_orientation,
        is.minmax_stackability
        )
end

function ItemizedStack(supplier_order, supplier_dock_order, plant_dock_order)
    return ItemizedStack(supplier_order, supplier_dock_order, plant_dock_order, :Free)
end


function ItemizedStack(os::OrderedStack)
    return ItemizedStack(os, :Free)
end

function ItemizedStack(os::OrderedStack, forced_orientation)
    return ItemizedStack(
        os,
        Item[],
        0.0,
        0.0,
        forced_orientation,
        Inf
        )
end

function ItemizedStack(supplier_order, supplier_dock_order, plant_dock_order, forced_orientation)
    return ItemizedStack(OrderedStack(supplier_order, supplier_dock_order, plant_dock_order), forced_orientation)
end

function set_ordered_stack!(is::ItemizedStack, os::OrderedStack)
    is.ordered_stack = os
    return
end

Base.hash(a::ItemizedStack, h::UInt) = hash(a.items, hash(:ItemizedStack, h))
Base.:(==)(a::ItemizedStack, b::ItemizedStack) = isequal(a.items, b.items) # TODO problem when 2 different stacks have an identical item
# Base.length(a::ItemizedStack) = length(a.items)



get_dim(is::ItemizedStack) = get_dim(is.ordered_stack)

function get_height(is::ItemizedStack)

    return is.height
end

function set_height(truck::Truck, h)
    return Truck(
        truck.dim,
        h,
        truck.max_stack_density,
        truck.max_stack_weight,
        truck.supplier_orders,
        truck.supplier_dock_orders,
        truck.plant_dock_orders,
    
        truck.CM,
        truck.CJ_fm,
        truck.CJ_fc,
        truck.CJ_fh,
        truck.EM,
        truck.EJ_hr,
        truck.EJ_cr,
        truck.EJ_eh,
        truck.EM_mr,
        truck.EM_mm,
    )
end

function get_weight(is::ItemizedStack)

    return is.weight
end

get_stackability_code(is::ItemizedStack) = isempty(is.items) ? nothing : get_stackability_code(is.items[1])
get_plant(is::ItemizedStack) =  isempty(is.items) ? nothing : get_plant(is.items[1])
get_plant_dock(is::ItemizedStack) = isempty(is.items) ? nothing : get_plant_dock(is.items[1])
get_supplier(is::ItemizedStack) = isempty(is.items) ? nothing : get_supplier(is.items[1])
get_supplier_dock(is::ItemizedStack) = isempty(is.items) ? nothing : get_supplier_dock(is.items[1])
get_ordered_stack(is::ItemizedStack) = isnothing(is.ordered_stack) ? nothing : is.ordered_stack
get_items(is::ItemizedStack) = is.items

get_supplier_order(is::ItemizedStack) = isnothing(is.ordered_stack) ? error("Can't get supplier order of nothing") : get_supplier_order(is.ordered_stack)
get_supplier_dock_order(is::ItemizedStack) = isnothing(is.ordered_stack) ? error("Can't get supplier dock order of nothing") : get_supplier_dock_order(is.ordered_stack)
get_plant_dock_order(is::ItemizedStack) = isnothing(is.ordered_stack) ? error("Can't get plant dock order of nothing") : get_plant_dock_order(is.ordered_stack)

get_orders(is::ItemizedStack) = get_orders(is.ordered_stack)
get_forced_orientation(is::ItemizedStack) = is.forced_orientation

get_pos(is::ItemizedStack) = isnothing(is.ordered_stack) ? nothing : get_pos(is.ordered_stack)
function get_dim(is::ItemizedStack) 
    if isnothing(is.ordered_stack) || isnothing(get_dim(is.ordered_stack))
        if !isempty(get_items(is))
            return get_dim(get_items(is)[1])
        else
            return nothing
        end
    else 
        return get_dim(is.ordered_stack)
    end
end

get_minmax_stackability(is::ItemizedStack) = is.minmax_stackability

Base.show(io::IO, is::ItemizedStack) = 
    print(io, "ItemizedStack(", readable(get_pos(is)), ", ", readable(get_dim(is)), ", orders=", get_orders(is), ", ", length(get_items(is)), " item(s), height=", round(get_height(is), digits=3), ", weight=", round(get_weight(is), digits=3), ", ", get_forced_orientation(is), ")")


function newStack(is::ItemizedStack, p::Pos, d::Dim)
    is.ordered_stack = newStack(is.ordered_stack, p, d)
    return is
end

"""
    add_item!(is::ItemizedStack, it::Item)

Add item `it` to stack `is`, updating weight, height, and maybe orientation.
"""
function add_item!(is::ItemizedStack, it::Item)
    push!(get_items(is), it)
    is.weight += get_weight(it)
    is.height += get_height(it)
    is.minmax_stackability = min(is.minmax_stackability, get_max_stackability(it))
    # update forced orientation
    if get_forced_orientation(it) != :Free
        if get_forced_orientation(it) != get_forced_orientation(is) && get_forced_orientation(is) != :Free
            error("Stack already has a different forced orientation")
        end
        is.forced_orientation = get_forced_orientation(it)
    end
    if length(get_items(is)) == 1
        is.height += get_nesting_height(it)
        # update stackability code
        # is.stackability_code = get_stackability_code(it) # no need
        # no need to update dims because easily computable
    end
end


"""
    can_be_placed(solution, o::Pos, s::ItemizedStack, W, orientation::Symbol; precision=3, verbose=false)

Return true if stack `s` can be placed in partial solution `solution` without breaking any constraint.
"""
function can_be_placed(solution, o::Pos, s::ItemizedStack, truck::Truck, orientation::Symbol; precision=3, verbose=false)

    # check weight constraints
    # TODO provide valid_axle_pressure the stack `s` with position pos
    # instantiate pos and dim
    os = set_stack(get_ordered_stack(s), Stack(o, get_dim(s)))
    set_ordered_stack!(s, os)
    weight_constraint = valid_axle_pressure(collect(values(solution)), s, truck; fastexit=true, precision=precision)

    return weight_constraint && can_be_placed(solution, o, get_ordered_stack(s), truck, orientation; precision=precision, verbose=verbose)
end

"""
    valid_stack(s, it, max_height)

Return true if item `it` can be added to stack `s` without breaking any dynamic constraint.
"""
function valid_stack(stacks, s, it, truck; fastexit=false, precision=3, verbose=false)

    # tm_t, ej_e, ej_r, em_h, em_r, em_m = dist_stacks_to_trailer(stacks, s, get_weight(it), truck)
    if verbose
        println(1, " ", leqtol(get_height(s) + get_height(it), get_height(truck), precision))
        println(2, " ", leqtol(get_weight(s) + get_weight(it), get_max_weight(it), precision))
        println(3, " ", length(get_items(s)) <= get_minmax_stackability(s))
        println(4, " ", get_forced_orientation(it) == :Free || get_forced_orientation(s) == :Free || get_forced_orientation(s) == get_forced_orientation(it))
        println(5, " ", leqtol((get_weight(s) + get_weight(it))/(get_dim(s).le * get_dim(s).wi), get_max_stack_density(truck), precision))
        println(6, " ", leqtol(get_weight(s) + get_weight(it), get_max_stack_weight(truck), precision))
        println(7, " ", valid_axle_pressure(stacks, s, it, truck; fastexit=fastexit, precision=precision))
    end
    return leqtol(get_height(s) + get_height(it), get_height(truck), precision) && 
    leqtol(get_weight(s) + get_weight(it), get_max_weight(it), precision) && # TODO might remove get_weight(it) since we want to limit the weight on bottom item
    length(get_items(s)) <= get_minmax_stackability(s) && # we need to find the smallest max_stackability of the pile
    (get_forced_orientation(it) == :Free || get_forced_orientation(s) == :Free || get_forced_orientation(s) == get_forced_orientation(it)) &&
    leqtol((get_weight(s) + get_weight(it))/(get_dim(s).le * get_dim(s).wi), get_max_stack_density(truck), precision) && # check density
    leqtol(get_weight(s) + get_weight(it), get_max_stack_weight(truck), precision) && # check max weight
    valid_axle_pressure(stacks, s, it, truck; fastexit=fastexit, precision=precision)
end

function valid_axle_pressure(stacks, stack::ItemizedStack, item::Item, truck; fastexit=false, precision=precision)

    dummystack = ItemizedStack(get_ordered_stack(stack)) # retrieve same position

    add_item!(dummystack, item)

    return valid_axle_pressure(stacks, dummystack, truck; fastexit=fastexit, precision=precision)
end


function valid_axle_pressure(stacks, s::ItemizedStack, truck::Truck; fastexit=false, precision=3, verbose=false)
    if verbose
        println("enter valid_axle_pressure")
    end
    # take stacks of all suppliers at first
    sortedsuppliers = sort(get_suppliers(truck), by= x -> get_supplier_orders(truck)[x])
    filteredstacks = copy(stacks)
    
    # display(sortedsuppliers)
    while !isempty(sortedsuppliers)

        if fastexit
            if !(get_supplier(s) in sortedsuppliers)
                break
            end
        end

        # find corresponding stacks
        filteredstacks = filter(x -> get_supplier(x) in sortedsuppliers, filteredstacks)
        
        # calculate pressure on middle and rear axles
        tm_t, ej_e, ej_r, em_h, em_r, em_m = nothing, nothing, nothing, nothing, nothing, nothing
        if get_supplier(s) in sortedsuppliers
            tm_t, ej_e, ej_r, em_h, em_r, em_m = dist_stacks_to_trailer(filteredstacks, s, truck)
        else
            tm_t, ej_e, ej_r, em_h, em_r, em_m = dist_stacks_to_trailer(filteredstacks, truck)
        end
        if verbose
            display(sortedsuppliers)
            display((tm_t, ej_e, ej_r, em_h, em_r, em_m))
        end

        if !(leqtol(em_m, get_EM_mm(truck), precision) && # max weight on middle axle
            leqtol(em_r, get_EM_mr(truck), precision)) # max weight on rear axle
            return false
        end
        # if ok, remove the last supplier
        pop!(sortedsuppliers)
        # do it again
    
    end
    return true
end


# function dist_stacks_to_trailer(allstacks, stack, item, truck)

#     dummystack = ItemizedStack(get_ordered_stack(stack)) # retreive same position

#     add_item!(dummystack, item)

#     # allstacks_copy = copy(allstacks)
#     # push!(allstacks_copy, dummystack)

#     return dist_stacks_to_trailer(allstacks, dummystack, truck)

# end

function dist_stacks_to_trailer(allstacks, stack, truck)


    allstacks_copy = copy(allstacks)
    push!(allstacks_copy, stack)

    return dist_stacks_to_trailer(allstacks_copy, truck)

end

function dist_stacks_to_trailer(allstacks, truck)

    tm_t = isempty(allstacks) ? 0.0 : sum(get_weight(s) for s in allstacks)

    ej_e = isempty(allstacks) ? 0 : sum((get_pos(s).x + get_dim(s).le/2) * get_weight(s) for s in allstacks) / tm_t

    ej_r = get_EJ_eh(truck) + get_EJ_hr(truck) - ej_e

    em_h = (tm_t * ej_r + get_EM(truck) * get_EJ_cr(truck)) / get_EJ_hr(truck)

    em_r = tm_t + get_EM(truck) - em_h

    em_m = (get_CM(truck) * get_CJ_fc(truck) + em_h * get_CJ_fh(truck)) / get_CJ_fm(truck)

    return tm_t, ej_e, ej_r, em_h, em_r, em_m

end


"""
    is_candidate_stack(stack, it)

Return true if stack `s` has the same loading order as `it` and the same stackability code.
"""
function is_candidate_stack(stack, it)
    return get_supplier(stack) == get_supplier(it) && 
    get_supplier_dock(stack) == get_supplier_dock(it) && 
    get_plant_dock(stack) == get_plant_dock(it) && 
    get_stackability_code(stack) == get_stackability_code(it)
end

"""
    find_candidate_stacks(it, stacks)

Return candidate stacks in `stacks` in regard to item `it`.
"""
function find_candidate_stacks(it, stacks)
    return filter(
        x -> is_candidate_stack(x, it), stacks)
end

"""
    make_stacks(items::Vector{Item}, plant_dock_orders, supplier_orders, supplier_dock_orders, max_height)

Return stacks made up of input `items` as to satisfy loading orders and stackability code aswell as dynamic constraints such
as height, weight, etc...

Items are supposed to have the same truck and thus the same plant.
"""
function make_stacks(items::Vector{Item}, plant_dock_orders, supplier_orders, supplier_dock_orders, max_height)
    # TODO check all items have same plant

    # sort item indices by supplier, supplier_dock, plant, plant_dock, stackability_code...
    # or not

    # build stacks online
    stacks = Dict{Any, Vector{ItemizedStack}}() # use supplier as key to group similar stacks TODO find better key
    for it in items
        if !(get_supplier(it) in keys(stacks))
            stacks[get_supplier(it)] = []
        end
        candidate_stacks = stacks[get_supplier(it)]
        
        candidate_stacks = find_candidate_stacks(it, candidate_stacks)

        found_stack = false
        # Choose a stack if height is ok and check max_weight too + max stackability
        # also check orientation
        for s in candidate_stacks
            if valid_stack(s, it, max_height)
                add_item!(s, it)
                found_stack = true
                break
            end
        end
        if !found_stack
            # create new stack
            # with good load orders
            newstack = ItemizedStack(
                supplier_orders[get_supplier(it)], 
                supplier_dock_orders[get_supplier(it)][get_supplier_dock(it)],
                plant_dock_orders[get_plant_dock(it)]
                )
            add_item!(newstack, it)

            push!(stacks[get_supplier(it)], newstack) 
        end
    end
    return stacks
end

"""
    rand_products(min_products, max_products, max_weight, max_items_per_stack)

Generate a random number between min_products and max_products of products of weight randomly taken between 0 and `max_weight` and
max items between 1 and `max_items_per_stack`.
"""
function rand_products(min_products, max_products, max_weight, max_items_per_stack)
    products = Vector{Product}(undef, rand(min_products:max_products))
    for i in 1:length(products)
        products[i] = Product(rand(1:max_items_per_stack), rand() * max_weight)
    end
    return products
end

"""
    rand_items(n, products, max_height, max_weight, L, W, plant; min_dim=0.001)

Generate `n` items of random properties. 

# Arguments
- `n`: the number of items to generate.
- `products`: products from which items should be generated.
- `max_height`: maximum height of an item.
- `max_weight`: maximum weight of an item.
- `L`: length of a truck. # TODO replace with maxlength
- `H`: height of a truck. # TODO same
- `plant`: code for the plant.
- `min_dim=0.001`: minimum value of length of width of an item.
"""
function rand_items(n, products, max_height, max_weight, L, W, plant; min_dim=0.001)

    # Product(
    # max_stackability::Integer,
    # max_weight::Real
    # )

    # Item(
    # time_window::@NamedTuple{earliest::Any, latest::Any},
    # dim::Dim,
    # # pos::Pos,
    # height::Real,
    # weight::Real,
    # stackability_code,
    # forced_orientation::Symbol,
    # plant,
    # plant_dock,
    # supplier,
    # supplier_dock,
    # inventory_cost,
    # nesting_height,
    # product::Product
    # )

    # Allocate table
    items = Vector{Union{Missing, Item}}(missing, n)
    n_filled = 0

    plant_docks = []
    supplier_docks = Dict() # key: supplier

    # while there are slots to fill
    while n_filled < n

        # draw random number of similar items
        nb = rand(1:n-n_filled)

        # affect random product
        p = rand(products)
    
        # generate random properties
        earliest = rand(0:365)
        
        time_window = (earliest=earliest, latest=rand(earliest:365))
        dim = Dim(max(min_dim, rand() * L), max(min_dim, rand() * W))
        # pos::Pos
        height = max([1, rand() * max_height]...)
        weight = max([1, rand() * max_weight]...)
        stackability_code = randstring(8)
        forced_orientation = rand([:Free, :Horizontal, :Vertical])

        plant_dock = nothing
        # choose existing or new dock
        if rand([true, false]) || isempty(plant_docks)
            plant_dock = randstring(8)
            push!(plant_docks, plant_dock)
        else
            plant_dock = rand(plant_docks)
        end
        if rand([true, false]) || isempty(keys(supplier_docks))
            supplier = randstring(8)
            supplier_docks[supplier] = []
        else
            supplier = rand(keys(supplier_docks))
        end
        supplier_dock = nothing
        if rand([true, false]) || isempty(supplier_docks[supplier])
            supplier_dock = randstring(8)
            push!(supplier_docks[supplier], supplier_dock)
        else
            supplier_dock = rand(supplier_docks[supplier])
        end

        inventory_cost = rand(1:10)
        nesting_height = rand([0, rand() * 10])

        for i in n_filled+1:n_filled + nb
            items[i] = Item(
                time_window,
                dim,
                height,
                weight,
                stackability_code,
                forced_orientation,
                plant,
                plant_dock,
                supplier,
                supplier_dock,
                inventory_cost,
                nesting_height,
                p)
        end

        n_filled += nb

    end


    # return table

    return items
end

"""
Assign ordered_stacks to ItemizedStacks
"""
function combine!(stacks::Dict{Integer, Stack}, ordered_stacks::Vector{Pair{Integer, OrderedStack}})
    for (i, os) in ordered_stacks
        set_ordered_stack!(stacks[i], os)
    end

    return
end

"""

Take `stack::Stack` and generate items to return ItemizedStacks. Is used when 
randomly generating an instance of itemizedstacks reusing the simple stack instance generator.
"""
function itemize!(truck::Truck, stacks::Dict{Integer, Stack})::Vector{Pair{Integer, ItemizedStack}}
    # give loading orders to stacks
    ordered_stacks, nbsupplier, nbsupplier_dock, nbplant_dock = order_instance(stacks)

    # create suppliers, supplier docks and plant docks based on order of stacks
    suppliers = Dict()
    supplier_docks = Dict()
    plant_docks = Dict()
    for (i, os) in ordered_stacks
        suppliers[get_supplier_order(os)] = randstring(8)
        supplier_docks[get_supplier_dock_order(os)] = randstring(8)
        plant_docks[get_plant_dock_order(os)] = randstring(8)
    end

    # TODO make sure drawn names aren't the same
    supplier_supplier_dock = Dict()

    set_supplier_orders!(truck, Dict(p[2] => p[1] for p in suppliers))
    # set_supplier_dock_orders!(truck, Dict(p[2] => p[1] for p in supplier_docks))
    set_plant_dock_orders!(truck, Dict(p[2] => p[1] for p in plant_docks))

    supplier_dock_orders = Dict(p[2] => p[1] for p in supplier_docks)

    # list of itemizedstacks
    istacks = []

    # create a bunch of products
    products = rand_products(1, length(stacks), 100.0, 100)
    stackability_codes = Dict()
    # for each stack
    for (i, os) in ordered_stacks

        # link supplier dock to supplier
        if !(suppliers[get_supplier_order(os)] in keys(supplier_supplier_dock))
            supplier_supplier_dock[suppliers[get_supplier_order(os)]] = []
        end
        push!(supplier_supplier_dock[suppliers[get_supplier_order(os)]], supplier_docks[get_supplier_dock_order(os)])

        # choose a random product
        p = rand(products)

        # choose a random number of items that satisfies product
        n = rand(1:get_max_stackability(p))
        
        # choose a random weight that satisfies product
        w = (1 + rand() * (get_max_weight(p) - 1)) / n

        # choose a random nesting height
        nesting_height = rand([0, rand() * 10])

        # determine height of items
        item_h = (get_height(truck) - nesting_height) / n

        # assign random orientation between free and orientation of stack
        stack_orient = get_dim(os).le > get_dim(os).wi ? :Horizontal : :Vertical
        # forced_orientation = rand([])

        # choose random time window
        earliest = rand(0:365)
        
        time_window = (earliest=earliest, latest=rand(earliest:365))

        # give dim of stack
        dim = get_dim(os)

        # choose a random stackability code (and store it)
        stackability_code = nothing
        if rand() > 0.8
            candidate_codes = filter(p -> 
                        (p[2].le == dim.le && p[2].wi == dim.wi) || 
                        (p[2].le == dim.wi && p[2].wi == dim.le), stackability_codes)
            if !isempty(candidate_codes)
                stackability_code = rand(candidate_codes)[1]
            end
        end
        if !isnothing(stackability_code)
            stackability_code = randstring(8)
            stackability_codes[stackability_code] = dim
        end
        

        # retieve plant (ignored), supplier and dock names from orders
        plant = ""
        supplier = suppliers[get_supplier_order(os)]
        supplier_dock = supplier_docks[get_supplier_dock_order(os)]
        plant_dock = plant_docks[get_plant_dock_order(os)]

        # choose random inventory cost
        inv_cost = 1
        items = []
        # create items (with random subset that have forced orientation)
        rand_subset_nb = rand(0:n)
        # create stack
        stack = ItemizedStack(os)
        wsum = 0
        for i in 1:n
            orientation = nothing
            if i <= rand_subset_nb
                # forced
                orientation = stack_orient
            else
                # free
                orientation = :Free
            end
            wsum += w
            item = Item(
                time_window,
                dim,
                item_h,
                w,
                stackability_code,
                orientation,
                plant,
                plant_dock,
                supplier,
                supplier_dock,
                inv_cost,
                nesting_height,
                p)
            add_item!(stack, copy(item))
        end
        push!(istacks, Pair(i, copy(stack)))
    end

    set_supplier_dock_orders!(truck, Dict(
                                        supplier => Dict(
                                                supplier_dock => supplier_dock_orders[supplier_dock]
                                                for supplier_dock in supplier_supplier_dock[supplier])
                                        for supplier in values(suppliers)))



    return istacks
end


function placeitem!(solution::Dict{T, S}, truck::Truck, item::Item, corners::Vector{<:AbstractPos}; precision=3, verbose=false) where {T <: Integer, S <: AbstractStack}
    torem = AbstractPos[]
    toadd = AbstractPos[]

    ind = isnothing(solution) || isempty(solution) ? 1 : max(collect(keys(solution))...) + 1
    # first find an existing stack which works for the item
    
    candidate_stacks = find_candidate_stacks(item, [p[2] for p in solution])

    found_stack = false
    # Choose a stack if height is ok and check max_weight too + max stackability
    # also check orientation # TODO improve complexity of the operation of doing the same calculation on all stacks
    # println("for s in candidate_stacks")
    for s in candidate_stacks
        # print(".")
        if valid_stack(collect(values(solution)), s, item, truck; fastexit=true, precision=precision)
            add_item!(s, item)
            found_stack = true
            break
        end
    end
    # println()

    # if no valid stack found create new stack containing the item and place item with standard placestack! function 
    if !found_stack
        # create new stack
        # with good load orders
        newstack = ItemizedStack(
            get_supplier_orders(truck)[get_supplier(item)], 
            get_supplier_dock_orders(truck)[get_supplier(item)][get_supplier_dock(item)],
            get_plant_dock_orders(truck)[get_plant_dock(item)]
            ) # created with no position

        add_item!(newstack, item)
        # display(get_items(newstack))
        torem, toadd = placestack!(solution, truck, ind, newstack, corners; precision=precision, verbose=verbose, loading_order=true) # DONE check weight constraints when adding a new stack
        # solution[ind] = copy(newstack)
    end

    return torem, toadd
end