using AutoHashEquals
using Random

include("placement.jl")
include("ordered_stacks.jl")

@auto_hash_equals struct Product
    max_stackability::Integer
    max_weight::Real
end

get_max_stackability(p::Product) = p.get_max_stackability
get_max_weight(p::Product) = p.get_max_weight

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


@auto_hash_equals struct ItemizedStack <: AbstractOrderedStack
    ordered_stack::OrderedStack
    items::Vector{Item}
    weight::Union{Nothing, Real}
    height::Union{Nothing, Real}

end

function ItemizedStack(supplier_order, supplier_dock_order, plant_dock_order)
    return ItemizedStack(
        OrderedStack(nothing, nothing, supplier_order, supplier_dock_order, plant_dock_order),
        Item[],
        nothing,
        nothing
        )
end


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

get_stackability_code(is::ItemizedStack) = isempty(items) ? missing : get_stackability_code(items[1])
get_plant(is::ItemizedStack) =  isempty(items) ? missing : get_plant(items[1])
get_plant_dock(is::ItemizedStack) = isempty(items) ? missing : get_plant_dock(items[1])
get_supplier(is::ItemizedStack) = isempty(items) ? missing : get_supplier(items[1])
get_supplier_dock(is::ItemizedStack) =isempty(items) ? missing : get_supplier_dock(items[1])
get_items(is::ItemizedStack) = is.items

get_pos(is::ItemizedStack) = get_pos(is.ordered_stack)
get_dim(is::ItemizedStack) = get_dim(is.ordered_stack)

function add_item!(is::ItemizedStack, it::Item)
    push!(get_items(is), it)
    is.weight += get_weight(it)
    is.height += get_height(it)
    # update forced orientation
    if get_orientation(it) != :Free
        if get_orientation(it) != get_orientation(is) && get_orientation(is) != :Free
            error("Stack already has a different forced orientation")
        end
        is.orientation = get_orientation(it)
    end
    if length(get_items(is)) == 1
        is.height += get_nesting_height(it)
        # update stackability code
        is.stackability_code = get_stackability_code(it)
        # upd dims
    end
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
        
        filter!(
            x -> get_supplier(x) == get_supplier(it) && 
            get_supplier_dock(x) == get_supplier_dock(it) && 
            get_plant_dock(x) == get_plant_dock(it) && 
            get_stackability_code(x) == get_stackability_code(it), candidate_stacks)

            found_stack = false
            # Choose a stack if height is ok and check max_weight too + max stackability
            # also check orientation
            for s in candidate_stacks
                if get_height(s) + get_height(it) <= max_height && 
                    get_weight(s) + get_weight(it) <= max_weight(it) && 
                    length(get_items(s)) <= get_max_stackability(it) &&
                    (get_orientation(s) == :Free || get_orientation(s) == get_orientation(it))
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
                    supplier_dock_orders[get_supplier_dock(it)],
                    plant_dock_orders[get_plant_dock(it)]
                    )
                add_item!(newstack, it)

                push!(stacks[get_supplier(it)], newstack) 
            end
    end

end

function rand_items(n, min_products, max_products, max_height, max_weight, max_items_per_stack, L, W, plant; min_dim=0.001)

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
    products = Vector{Product}(undef, rand(min_products:max_products))
    for i in 1:length(products)
        products[i] = Product(rand(1:max_items_per_stack), rand() * max_weight)
    end

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