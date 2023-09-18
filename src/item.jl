using AutoHashEquals
using Random

include("placement.jl")
include("ordered_stacks.jl")

@auto_hash_equals struct Product
    max_stackability::Integer
    max_weight::Real
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
        it.product # no deep copy necessary
        )
end

mutable struct ItemizedStack <: AbstractOrderedStack
    ordered_stack::Union{Nothing, OrderedStack}
    items::Vector{Item}
    weight::Real
    height::Real
    forced_orientation::Symbol
    function ItemizedStack(ordered_stack::OrderedStack,
        items::Vector{Item},
        weight::Real,
        height::Real,
        forced_orientation::Symbol)
        if !(forced_orientation in [:Free, :Horizontal, :Vertical])
            throw(ArgumentError("forced_orientation should be taken among [:Free, :Horizontal, :Vertical].\nGot $forced_orientation."))
        end
        new(ordered_stack, items, weight, height, forced_orientation)
    end

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
        forced_orientation
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
Base.:(==)(a::ItemizedStack, b::ItemizedStack) = isequal(a.items, b.items)
# Base.length(a::ItemizedStack) = length(a.items)

get_dim(is::ItemizedStack) = get_dim(is.ordered_stack)

function get_height(is::ItemizedStack)
    if isnothing(is.height)
        height = sum(get_height(i) for i in is.items) - get_nesting_height(is.items[1])
        set_height!(is, height)
        return height
    else
        return is.height
    end
end

function get_weight(is::ItemizedStack)
    if isnothing(is.weight)
        weight = sum(get_weight(i) for i in is.items)
        set_weight!(is, weight)
        return weight
    else
        return is.weight
    end
end

get_stackability_code(is::ItemizedStack) = isempty(is.items) ? nothing : get_stackability_code(is.items[1])
get_plant(is::ItemizedStack) =  isempty(is.items) ? nothing : get_plant(is.items[1])
get_plant_dock(is::ItemizedStack) = isempty(is.items) ? nothing : get_plant_dock(is.items[1])
get_supplier(is::ItemizedStack) = isempty(is.items) ? nothing : get_supplier(is.items[1])
get_supplier_dock(is::ItemizedStack) =isempty(is.items) ? nothing : get_supplier_dock(is.items[1])
get_items(is::ItemizedStack) = is.items

get_supplier_order(is::ItemizedStack) = isnothing(is.ordered_stack) ? error("Can't get supplier order of nothing") : get_supplier_order(is.ordered_stack)
get_supplier_dock_order(is::ItemizedStack) = isnothing(is.ordered_stack) ? error("Can't get supplier dock order of nothing") : get_supplier_dock_order(is.ordered_stack)
get_plant_dock_order(is::ItemizedStack) = isnothing(is.ordered_stack) ? error("Can't get plant dock order of nothing") : get_plant_dock_order(is.ordered_stack)

get_orders(is::ItemizedStack) = get_orders(is.ordered_stack)
get_forced_orientation(is::ItemizedStack) = is.forced_orientation

get_pos(is::ItemizedStack) = isnothing(is.ordered_stack) ? nothing : get_pos(is.ordered_stack)
get_dim(is::ItemizedStack) = isnothing(is.ordered_stack) ? nothing : get_dim(is.ordered_stack)

Base.show(io::IO, is::ItemizedStack) = 
    print(io, "ItemizedStack(Pos(", round(get_pos(is).x, digits=3), ", ", round(get_pos(is).y, digits=3), "), Dim(", round(get_dim(is).le, digits=3), ", ", round(get_dim(is).wi, digits=3), "), orders=", get_orders(is), ", ", length(get_items(is)), " item(s), height=", round(get_height(is), digits=3), ", weight=", round(get_weight(is), digits=3), ", ", get_forced_orientation(is), ")")

function add_item!(is::ItemizedStack, it::Item)
    push!(get_items(is), it)
    is.weight += get_weight(it)
    is.height += get_height(it)
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
        # upd dims
    end
end

function valid_stack(s, it, max_height)
    return get_height(s) + get_height(it) <= max_height && 
    get_weight(s) + get_weight(it) <= get_max_weight(it) && 
    length(get_items(s)) <= get_max_stackability(it) &&
    (get_forced_orientation(it) == :Free || get_forced_orientation(s) == :Free || get_forced_orientation(s) == get_forced_orientation(it))
end

function is_candidate_stack(stack, it)
    return get_supplier(stack) == get_supplier(it) && 
    get_supplier_dock(stack) == get_supplier_dock(it) && 
    get_plant_dock(stack) == get_plant_dock(it) && 
    get_stackability_code(stack) == get_stackability_code(it)
end

function find_candidate_stacks(it, stacks)
    return filter(
        x -> is_candidate_stack(x, it), stacks)
end

"""
For a same truck
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

function rand_products(min_products, max_products, max_weight, max_items_per_stack)
    products = Vector{Product}(undef, rand(min_products:max_products))
    for i in 1:length(products)
        products[i] = Product(rand(1:max_items_per_stack), rand() * max_weight * max_items_per_stack)
    end
    return products
end

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

    # create a random number of products
    # products = Vector{Product}(undef, rand(min_products:max_products))
    # for i in 1:length(products)
    #     products[i] = Product(rand(1:max_items_per_stack), rand() * max_weight * max_items_per_stack)
    # end

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


function itemize(stacks::Dict{Integer, Stack}, H)::Vector{Pair{Integer, ItemizedStack}}
    ordered_stacks = order_instance(stacks)

    istacks = []

    # create a bunch of products
    products = rand_products(1, length(stacks), 100.0, 100)
    stackability_codes = Dict()
    # for each stack
    for (i, os) in ordered_stacks
        # choose a random product
        p = rand(products)

        # choose a random number of items that satisfies product
        n = rand(1:get_max_stackability(p))
        
        # choose a random weight that satisfies product
        w = (1 + rand() * (get_max_weight(p) - 1)) / n

        # choose a random nesting height
        nesting_height = rand([0, rand() * 10])

        # determine height of items
        item_h = (H - nesting_height) / n

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
        # ignore plant, supplier and dock names
        # choose random inventory cost
        inv_cost = 1
        items = []
        # create items (with random subset that have forced orientation)
        rand_subset_nb = rand(0:n)
        # create stack
        stack = ItemizedStack(os)
        for i in 1:n
            orientation = nothing
            if i <= rand_subset_nb
                # forced
                orientation = stack_orient
            else
                # free
                orientation = :Free
            end
            item = Item(
                time_window,
                dim,
                item_h,
                w,
                stackability_code,
                orientation,
                "",
                "",
                "",
                "",
                inv_cost,
                nesting_height,
                p)
            
            add_item!(stack, item)
        end
        push!(istacks, Pair(i, stack))
    end
    return istacks
end