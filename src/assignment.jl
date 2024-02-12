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

function solve_tsi(t_trucks, i_items, TR)
    # i in i_items is line index of TR
    # t in t_trucks is column index of TR
    # extra trucks are multiples of original planned trucks in t_trucks
    nbuniqitems = length(Set(get_id(i_item[2]) for i_item in i_items))
    item_index = Dict(get_id(item) => ((i-1) % nbuniqitems)+1 for (i, item) in i_items)


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
        # display([iit[1] for iit in i_items])
        display([length(item_dispatch[k]) for k in keys(item_dispatch)])
        display([iit[1] for k in keys(item_dispatch) for iit in item_dispatch[k]])
        println("Sum of items = $(length(i_items) + sum([length(item_dispatch[k]) for k in keys(item_dispatch)]))")
        readline()
        error("Number of items increases")
        ## New truck phase
        i_init, item_init = i_items[1]
        # display([TR[((t-1) % ntrucks)+1, i_init] for (t, truck) in t_trucks])
        # display([can_take(truck, [], item_init) for (t, truck) in t_trucks])
        # readline()
        # take cheapest truck that can take first item (with most constraints)
        candi_t, candi_truck = t_trucks[findfirst(
            # x -> TR[i_init, x[1]] && can_take(x[2], item_dispatch[i_init], item_init), t_trucks
            t_truck -> 
                TR[((t_truck[1]-1) % ntrucks)+1, i_init] && 
                can_take(t_truck[2], [], item_init), 
            t_trucks
        )]
        # remove the truck candidate from candidate trucks
        filter!(t_truck -> t_truck[1] != candi_t, t_trucks)
        
        # Add extra truck to candidate trucks
        push!(t_trucks, Pair(candi_t + ntrucks, candi_truck)) # TODO keep it sorted
        sort!(t_trucks, by=x -> (get_cost(x[2]), ((.*).(versatility(x[2]), -1))...))
        if !issorted(t_trucks, by=t_truck -> (get_cost(t_truck[2]), ((.*).(versatility(t_truck[2]), -1))...))
            display([get_id(t_truck[2]) for t_truck in t_trucks])
            display(map(t_truck -> (get_cost(t_truck[2]), ((.*).(versatility(t_truck[2]), -1))...), t_trucks))
            error("t_trucks is not sorted.")
        end

        item_dispatch[candi_t] = []
        
        ## Assigning to trucks phase
        # then for each item remaining to place
        for (c, (i, item)) in enumerate(i_items)
            # fill most versatile trucks first
            # TODO first fill new truck?
            for (t, truck) in 
                [p for truck_list in [[(candi_t, candi_truck)] , used_trucks] for p in truck_list]
                # if max weight is ok for considered truck and volume is still ok aswell
                if TR[((t-1) % ntrucks)+1, ((i-1) % nbuniqitems)+1] && can_take(truck, [x[2] for x in item_dispatch[t]], item_init)
                    # add item to list to add to the truck
                    push!(item_dispatch[t], Pair(i, item)) # excess items in i_items
                    break
                end
            end

            display_progress(c, length(i_items); name="Item assignment")
        end
        # remove asssigned items from i_items
        filter!(x -> !(x in [
            s for k in keys(item_dispatch) for s in item_dispatch[k]
            ]), i_items) # balance out assigned items
        notplaced_global = []
        push!(used_trucks, Pair(candi_t, candi_truck)) # TODO keep used_trucks sorted


        # solve subproblems, retrieve items which couldn't be placed
        for (c, (t, truck)) in enumerate(used_trucks)
            solution[t], notplaced = BLtruck([x[2] for x in item_dispatch[t]], truck)
            # remove items of stacks which poke out of truck
            append!(
                notplaced, 
                [
                    item 
                    for (s, stack) in solution[t] 
                    for item in get_items(stack)
                    if get_pos(stack).x + get_dim(stack).le > get_dim(truck).le
                ]
            )
            # clean up solution
            filter!((s_stack) -> 
                get_pos(s_stack[2]).x + get_dim(s_stack[2]).le <= get_dim(truck).le, solution[t]
            )
                
            # transfer notplaced to notplaced_global
            append!(notplaced_global, map(x -> Pair(item_index[get_id(x)], x), notplaced))
            
            # remove notplaced items from item_dispatch of the truck
            filter!(x -> !(x in notplaced), item_dispatch[t])
            display_progress(c, length(used_trucks); name="Solving trucks")
        end
        



        # Already done by design?
        if !issorted(used_trucks, by= t_truck -> 
            (get_cost(t_truck[2]), (.*).(versatility(t_truck[2]), -1)...))

            error("used_trucks is not sorted.")
        end
        
        append!(i_items, notplaced_global)
        sort!(i_items, by= i_item -> (.*).(constraints(i_item[2], ((i_item[1]-1) % nbuniqitems)+1, TR), -1)) # TODO insert in placed?
        
    end
    
    return used_trucks, solution
end