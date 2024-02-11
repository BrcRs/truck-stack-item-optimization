include("dim.jl")
include("truck.jl")
include("item.jl")

function versatility(truck)
    volume = get_dim(truck).le * get_dim(truck).wi * get_height(truck)
    return volume, get_TMm(truck), get_max_stack_density(truck), get_EM_mm(truck), get_EM_mr(truck)
end

function constraints(item, i, TR)
    f_orient = get_orientation(item) != :none ? 1 : 0
    trucks_available = sum(TR[:, i])
    le_wi_ratio = get_dim(item).le / get_dim(item).wi
    shape_weirdness = max(le_wi_ratio, inv(le_wi_ratio))
    volume = get_dim(item).le * get_dim(item).wi * get_height(item)
    width_time_window = get_time_window(item)[2] - get_time_window(item)[1]

    return trucks_available, width_time_window, volume, shape_weirdness, f_orient 
end

function can_take(truck, solution, item)
    # error("TODO check if same plant (just check TR basically?)")
    volume_left = get_volume(truck) - get_volume(solution)
    return get_volume(item) <= volume_left && 
        get_weight(solution) + get_weight(item) <= get_TMm(truck)
end

# TODO treat case when item is constantly candidate to a truck but weight or other 
# constraints makes itsd integration impossible

function solve_tsi(t_trucks, i_items, TR)
    # i in i_items is line index of TR
    # t in t_trucks is column index of TR
    # extra trucks are multiples of original planned trucks in t_trucks

    ntrucks = length(t_trucks)

    used_trucks = []

    solution = Dict{Any, Any}() # index is used trucks

    # sort trucks by increasing cost and by decreasing versatility
    sort!(t_trucks, by=x -> (get_cost(x[2]), -versatility(x[2])...))

    # sort items by decreasing constraints
    sort!(i_items, by= x -> -constraints(x[2], x[1], TR))

    # while there are still items to place
    while !isempty(i_items)
        ## New truck phase
        i_init, item_init = i_items[1]
        
        # take cheapest truck that can take first item (with most constraints)
        candi_t, candi_truck = findfirst(
            # x -> TR[i_init, x[1]] && can_take(x[2], solution[i_init], item_init), t_trucks
            (t, truck) -> TR[((t-1) % ntrucks)+1, i_init] && can_take(truck, [], item_init), t_trucks
        )
        
        # Add extra truck to candidate trucks
        push!(t_trucks, Pair(candi_t + ntrucks, candi_truck)) # TODO keep it sorted
        if !issorted(t_trucks, by=t, truck -> (get_cost(truck), -versatility(truck)...))
            error("t_trucks is not sorted.")
        end
        
        
        ## Assigning to trucks phase
        # then for each item remaining to place
        for (i, item) in i_items
            # fill most versatile trucks first
            # TODO first fill new truck?
            for (t, truck) in 
                [p for p in truck_list for truck_list in [[(candi_t, candi_truck)], used_trucks]]
                # if max weight is ok for considered truck and volume is still ok aswell
                if TR[((t-1) % ntrucks)+1, i] && can_take(truck, solution[t], item_init)
                    # add item to list to add to the truck
                    push!(solution[t], item)
                end
            end
        end
        filter!(x -> !(x[2] in [s for k in keys(solution) for s in solution[k]]), i_items)
        notplaced_global = []
        # solve subproblems, retrieve items which couldn't be placed
        for (t, truck) in used_trucks
            solution, notplaced = BLtruck(solution[t], truck)
            append!(notplaced_global, notplaced)
        end
        
        push!(used_trucks, Pair(candi_t, candi_truck)) # TODO keep used_trucks sorted
        # Already done by design?
        if !issorted(used_trucks, by= t, truck -> (get_cost(truck), -versatility(truck)...))
            error("used_trucks is not sorted.")
        end
        
        append!(i_items, notplaced_global)
        sort!(i_items, by= i, item -> -constraints(item, i, TR)) # TODO insert in placed?
    end
    
    return used_trucks, solution
end