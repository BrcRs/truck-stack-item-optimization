include("dim.jl")
include("truck.jl")
include("item.jl")
include("progress.jl")
include("placement_algorithms.jl")

function versatility(truck)
    volume = get_dim(truck).le * get_dim(truck).wi * get_height(truck)
    return volume, get_TMm(truck), get_max_stack_density(truck), get_EM_mm(truck), get_EM_mr(truck)
end

function constraints(item, i, TR)
    f_orient = get_forced_orientation(item) != :none ? 1 : 0
    trucks_available = sum(TR[:, i])
    le_wi_ratio = get_dim(item).le / get_dim(item).wi
    shape_weirdness = max(le_wi_ratio, inv(le_wi_ratio))
    volume = get_dim(item).le * get_dim(item).wi * get_height(item)
    width_time_window = get_time_window(item)[2] - get_time_window(item)[1]

    return trucks_available, width_time_window, volume, shape_weirdness, f_orient 
end

function can_take(truck, items, it)
    # error("TODO check if same plant (just check TR basically?)")
    volume_left = get_volume(truck) - (isempty(items) ? 0 : sum(get_volume(s) for s in items))
    return get_volume(it) <= volume_left && 
    (isempty(items) ? 0 : sum(get_weight(s) for s in items)) + get_weight(it) <= get_TMm(truck)
end

# TODO treat case when item is constantly candidate to a truck but weight or other 
# constraints makes itsd integration impossible

function assigning_to_truck_step!(i, item, t, truck, TR, item_dispatch, ntrucks, nbuniqitems)
    # fill most versatile trucks first
    # TODO first fill new truck?
    # if max weight is ok for considered truck and volume is still ok aswell
    if TR[((t-1) % ntrucks)+1, ((i-1) % nbuniqitems)+1] && can_take(truck, [x[2] for x in item_dispatch[t]], item)
        # add item to list to add to the truck
        push!(item_dispatch[t], Pair(i, item)) # excess items in i_items
        return true
    end
    return false
end

truck_sort_fn(t_truck) = (get_cost(t_truck[2]), ((.*).(versatility(t_truck[2]), -1))...)

item_sort_fn(i_item, nbuniqitems, TR) = (.*).(constraints(i_item[2], ((i_item[1]-1) % nbuniqitems)+1, TR), -1)

function valid_truck(t_truck, i_item, ntrucks, nbuniqitems, solution_t, TR)

    t_mod = ((t_truck[1]-1) % ntrucks)+1
    i_mod = ((i_item[1]-1) % nbuniqitems)+1

    return TR[t_mod, i_mod] && can_take(t_truck[2], solution_t, i_item[2])
end

function select_new_truck!(t_trucks, i_item, ntrucks, nbuniqitems, TR)
    ## New truck phase
    # take cheapest truck that can take first item (with most constraints)
    candi_t, candi_truck = t_trucks[findfirst(
        t_truck -> valid_truck(t_truck, i_item, ntrucks, nbuniqitems, [], TR),
        t_trucks
    )]
    # remove the truck candidate from candidate trucks
    filter!(t_truck -> t_truck[1] != candi_t, t_trucks)
    
    
    # Add extra truck to candidate trucks
    push!(t_trucks, Pair(candi_t + ntrucks, candi_truck)) # TODO keep it sorted
    sort!(t_trucks, by=truck_sort_fn)
    if !issorted(t_trucks, by=truck_sort_fn)
        display([get_id(t_truck[2]) for t_truck in t_trucks])
        display(map(truck_sort_fn, t_trucks))
        error("t_trucks is not sorted.")
    end
    

    return candi_t, candi_truck

end

function assign_to_trucks!(i_items, candi_t_truck, used_trucks, TR, item_dispatch, ntrucks, nbuniqitems)

    ## Assigning to trucks phase
    # then for each item remaining to place
    for (c, (i, item)) in enumerate(i_items)
        tmp_trucklist = [p for truck_list in [[candi_t_truck], used_trucks] for p in truck_list]
        for (t, truck) in tmp_trucklist
            if assigning_to_truck_step!(i, item, t, truck, TR, item_dispatch, ntrucks, nbuniqitems)
                break
            end
        end
        display_progress(c, length(i_items); name="Item assignment")
    end

    # remove asssigned items from i_items
    flat_item_dispatch = [s for k in keys(item_dispatch) for s in item_dispatch[k]]
    filter!(x -> !(x in flat_item_dispatch), i_items) # balance out assigned items
    push!(used_trucks, candi_t_truck) # TODO keep used_trucks sorted
end

function solve_tsi_step!(item_dispatch, used_trucks, solution, item_index)
    
    notplaced_global = Pair{Int64, Item}[]
    # solve subproblems, retrieve items which couldn't be placed
    for (c, (t, truck)) in enumerate(used_trucks)
        count_before = length(item_dispatch[t])
        solution[t], notplaced = BLtruck([x[2] for x in item_dispatch[t]], truck)
        # println("item dispatch")
        # println([string(get_id(x[2]), "::", get_copy_number(x[2])) for x in item_dispatch[t]])
        # readline() 
        ########### DEBUG
        # remove items of stacks which poke out of truck
        # println("not placed")
        # println([string(get_id(x), "::", get_copy_number(x)) for x in notplaced])
        # readline()

        append!(
            notplaced, 
            [
                item 
                for (s, stack) in solution[t] 
                for item in get_items(stack)
                if get_pos(stack).x + get_dim(stack).le > get_dim(truck).le # when commented no more pb?
            ]
        )
        # println("not placed extended")
        # println([string(get_id(x), "::", get_copy_number(x)) for x in notplaced])
        # readline()
        # the bug does not appear when notplaced is either not altered
        # or contains all items
        ###################### This block deletes items on accident

        # clean up solution
        filter!((s_stack) -> 
            get_pos(s_stack[2]).x + get_dim(s_stack[2]).le <= get_dim(truck).le, solution[t]
        )
        
        # transfer notplaced to notplaced_global
        append!(notplaced_global, map(item -> Pair(item_index[string(get_id(item), "::", get_copy_number(item))], item), notplaced))
        # println("not placed global")
        # println([string(get_id(x[2]), "::", get_copy_number(x[2])) for x in notplaced_global])
        # readline()
        # remove notplaced items from item_dispatch of the truck
        filter!(x -> !(string(get_id(x[2]), "::", get_copy_number(x[2])) in [string(get_id(y), "::", get_copy_number(y)) for y in notplaced]), item_dispatch[t])
        # println("filtered item_dispatch")
        # println([string(get_id(x[2]), "::", get_copy_number(x[2])) for x in item_dispatch[t]])
        # readline()


        count_after = length(item_dispatch[t]) + length(notplaced)
        if count_before != count_after
            error("Items disappeared/duplicated\n$count_before != $count_after")
        end
        display_progress(c, length(used_trucks); name="Solving trucks")
    end
    

    # Already done by design?
    if !issorted(used_trucks, by=truck_sort_fn)

        error("used_trucks is not sorted.")
    end
    
    return notplaced_global

end

function solve_tsi(t_trucks, i_items, TR)
    # i in i_items is line index of TR
    # t in t_trucks is column index of TR
    # extra trucks are multiples of original planned trucks in t_trucks

    nbuniqitems = length(Set(get_id(i_item[2]) for i_item in i_items))
    item_index = Dict(string(get_id(item), "::", get_copy_number(item)) => i for (i, item) in i_items)


    ntrucks = length(t_trucks)

    used_trucks = []

    item_dispatch = Dict{Any, Vector{Pair{Integer, Item}}}() # index is used trucks
    solution = Dict()
    # sort trucks by increasing cost and by decreasing versatility
    sort!(t_trucks, by=x -> (get_cost(x[2]), ((.*).(versatility(x[2]), -1))...))

    # sort items by decreasing constraints
    # ((x[1]-1) % nbuniqitems)+1
    sort!(i_items, by= x -> ((.*).(constraints(x[2], ((x[1]-1) % nbuniqitems)+1, TR), -1)))

    # while there are still items to place
    while !isempty(i_items)
        println("Items left: $(length(i_items))")
        # println([iit[1] for iit in i_items])
        # display([length(item_dispatch[k]) for k in keys(item_dispatch)])
        # # display([iit[1] for k in keys(item_dispatch) for iit in item_dispatch[k]])
        # println("Sum of items = $(length(i_items) + sum([length(item_dispatch[k]) for k in keys(item_dispatch)]))")
        # readline()
        # error("Number of items increases")
        candi_t, candi_truck = select_new_truck!(t_trucks, i_items[1], ntrucks, nbuniqitems, TR)

        item_dispatch[candi_t] = []
        
        # assign to truck shouldn't make duplicates
        assign_to_trucks!(i_items, Pair(candi_t, candi_truck), used_trucks, TR, item_dispatch, ntrucks, nbuniqitems)
        
        notplaced_global = solve_tsi_step!(item_dispatch, used_trucks, solution, item_index)
        for stuff in notplaced_global
            if stuff in i_items
                error("$stuff in i_items")
            end
        end

        append!(i_items, notplaced_global)

        sort!(i_items, by= i_item -> item_sort_fn(i_item, nbuniqitems, TR)) # TODO insert in placed?
    end
    
    return used_trucks, solution
end