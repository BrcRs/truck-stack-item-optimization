using JuMP
using LinearAlgebra
using SparseArrays 
using Dates

include("dim.jl")
include("truck.jl")
include("item.jl")
include("progress.jl")
include("placement_algorithms.jl")
include("placement_visualizer.jl")

function versatility(truck)
    volume = get_dim(truck).le * get_dim(truck).wi * get_height(truck)
    return -volume, -get_TMm(truck), -get_max_stack_density(truck), -get_EM_mm(truck), -get_EM_mr(truck)
end

function constraints(item, i, TR)
    f_orient = get_forced_orientation(item) != :none ? 1 : 0
    trucks_available = sum(TR[:, i])
    le_wi_ratio = get_dim(item).le / get_dim(item).wi
    shape_weirdness = max(le_wi_ratio, inv(le_wi_ratio))
    volume = get_dim(item).le * get_dim(item).wi * get_height(item)
    width_time_window = get_time_window(item)[2] - get_time_window(item)[1]

    # return trucks_available, width_time_window, -volume, shape_weirdness, -f_orient 
    return get_plant(item), get_plant_dock(item), get_stackability_code(item), get_time_window(item)[1], get_time_window(item)[2],
        get_supplier(item), get_supplier_dock(item), -get_dim(item).le, -get_dim(item).wi, -get_height(item)
end

function can_take(truck, it, items_volume, items_weight)
    # error("TODO check if same plant (just check TR basically?)")
    # volume_left = get_volume(truck) - (isempty(items) ? 0 : sum(get_volume(s) for s in items))
    volume_left = get_volume(truck) - items_volume
    return get_volume(it) <= volume_left && 
    # (isempty(items) ? 0 : sum(get_weight(s) for s in items)) + get_weight(it) <= get_TMm(truck)
    items_weight + get_weight(it) <= get_TMm(truck)
end

# TODO treat case when item is constantly candidate to a truck but weight or other 
# constraints makes its integration impossible

function assigning_to_truck_step!(ind, i, item, t, truck, TR, item_dispatch, ntrucks, 
    nbuniqitems, items_volume, items_weight)
    # fill most versatile trucks first
    # TODO first fill new truck?
    # if max weight is ok for considered truck and volume is still ok aswell
    # if TR[((t-1) % ntrucks)+1, ((i-1) % nbuniqitems)+1] && can_take(truck, [x[2] for x in item_dispatch[t]], item)
    if valid_truck(Pair(t, truck), Pair(i, item), ntrucks, nbuniqitems, items_volume[t], items_weight[t], TR)
        # add item to list to add to the truck
        # push!(item_dispatch[t], Pair(i, item))
        item_dispatch[ind] = t
        items_volume[t] += get_volume(item)
        items_weight[t] += get_weight(item)
        return true
    end
    return false
end
#TODO multiply sum TR by numbers of copies of each item
function truck_sort_fn(t_truck, ntrucks, item_dispatch, TR, compatible_items)
    # nbcompatible_items = count(x -> item_dispatch[x] == -1, compatible_items[t_mod(t_truck[1], ntrucks)])
    nbcompatible_items = length(compatible_items[t_mod(t_truck[1], ntrucks)])

    return -nbcompatible_items, t_truck[1] <= ntrucks ? 0. : get_cost(t_truck[2]), -sum(TR[((t_truck[1]-1) % ntrucks)+1, :]), versatility(t_truck[2])...
end
t_mod(t, ntrucks) = ((t-1) % ntrucks)+1
i_mod(i, nbuniqitems) = ((i-1) % nbuniqitems)+1

item_sort_fn(i_item, nbuniqitems, TR) = constraints(i_item[2], ((i_item[1]-1) % nbuniqitems)+1, TR)

function valid_truck(t_truck, i_item, ntrucks, nbuniqitems, items_volume, items_weight, TR)

    t_mod = ((t_truck[1]-1) % ntrucks)+1
    i_mod = ((i_item[1]-1) % nbuniqitems)+1

    return TR[t_mod, i_mod] && can_take(t_truck[2], i_item[2], items_volume, items_weight)
end

function new_truck_id(t, truck, ntrucks)
    truck_cp_number = div(t-1, ntrucks)
    if truck_cp_number == 0
        return string("P", get_id(truck)[2:end])
    else
        return string("Q", split(get_id(truck)[2:end], "_")[1], "_", truck_cp_number)
    end
end

function select_new_truck!(t_trucks, i_item, ntrucks, nbuniqitems, TR, compatible_items, item_dispatch)
    ## New truck phase
    # take cheapest truck that can take first item (with most constraints)
    candi_t, candi_truck = t_trucks[findfirst(
        t_truck -> valid_truck(t_truck, i_item, ntrucks, nbuniqitems, 0., 0., TR),
        t_trucks
    )]
    # remove the truck candidate from candidate trucks
    # filter!(t_truck -> t_truck[1] != candi_t, t_trucks) #TODO maybe too long?
    deleteat!(t_trucks, findfirst(t_truck -> t_truck[1] == candi_t, t_trucks)) #TODO do at the end of select_multiple_trucks?
    # candi_truck = copy(candi_truck)
    candi_truck = set_id(candi_truck, new_truck_id(candi_t, candi_truck, ntrucks)) # change P to Q if extra + _ + copy number
    
    # Add extra truck to candidate trucks
    extra_t_truck = Pair(candi_t + ntrucks, candi_truck)
    # push!(t_trucks, extra_t_truck) # DONE keep it sorted
    insert_ind = findfirst(t_truck -> truck_sort_fn(extra_t_truck, ntrucks, item_dispatch, TR, compatible_items) < truck_sort_fn(t_truck, ntrucks, item_dispatch, TR, compatible_items), t_trucks)
    insert!(
        t_trucks, 
        isnothing(insert_ind) ? length(t_trucks) + 1 : insert_ind, 
        extra_t_truck
    )
    # push!(t_trucks, extra_t_truck)
    # partialsort!(
    #     t_trucks, 
    #     length(t_trucks); 
    #     by=t_truck -> truck_sort_fn(extra_t_truck, ntrucks, TR)
    # )
    

    return candi_t, candi_truck

end

function assign_to_trucks!(
    i_items, items_indices, candi_t_truck_list, used_trucks, TR, item_dispatch, ntrucks, 
    nbuniqitems, compatible_trucks, compatible_items, item_undispatch; limit=nothing
)

    tmp_trucklist = [p for truck_list in [candi_t_truck_list, used_trucks] for p in truck_list]
    # start = time()
    items_volume = Dict(
        t => isnothing(findfirst(x -> item_dispatch[x] == t, 1:length(i_items))) ?
            0. :
            sum(
                [
                    get_volume(i_items[ind][2]) 
                    for ind in filter(x -> item_dispatch[x] == t, 1:length(i_items))
                ]
            ) 
        for (t, truck) in tmp_trucklist
    )
    items_weight = Dict(
        t => isnothing(findfirst(x -> item_dispatch[x] == t, 1:length(i_items))) ?
            0. :
            sum(
                [
                    get_weight(i_items[ind][2]) 
                    for ind in filter(x -> item_dispatch[x] == t, 1:length(i_items))
                ]
            ) 
        for (t, truck) in tmp_trucklist 
    )
    # println("Time to compute weights and volumes: $(time() - start)")
    # DONE function too long to process
    # DONE Idea: keep in memory a list of truck indices for each "group of items" for which those trucks work
    # Then, only probe those trucks for the items instead of trying everything
    # DONE Keep count of trucks which rejected item i. Do not try to fit item i in this truck again 
    ## Assigning to trucks phase
    # then for each item remaining to place
    println()
    for (c, ind) in enumerate(items_indices)
        i, item = i_items[ind]
        if !isnothing(limit) && c > limit
            break
        end
        for (t, truck) in tmp_trucklist
            if !(((t-1) % ntrucks)+1 in compatible_trucks[((i-1) % nbuniqitems)+1]) ||
                t in item_undispatch[ind]
                continue
            end
            if assigning_to_truck_step!(ind, i, item, t, truck, TR, item_dispatch, ntrucks, nbuniqitems,
                items_volume, items_weight)
                break
            end
        end
        display_progress(
            c, 
            isnothing(limit) ? length(items_indices) : min(length(items_indices), limit); 
            name="Item assignment", 
            rate=5
        )
    end

    # remove asssigned items from i_items
    # flat_item_dispatch = [s for k in keys(item_dispatch) for s in item_dispatch[k]]
    # filter!(x -> !(x in flat_item_dispatch), i_items) # balance out assigned items
    
    append!(used_trucks, candi_t_truck_list)
    sort!(used_trucks, by=t_truck -> truck_sort_fn(t_truck, ntrucks, item_dispatch, TR, compatible_items))
    
    # append!(used_trucks, candi_t_truck_list)
    # partialsort!(
    #     used_trucks, 
    #     (length(used_trucks) - length(candi_t_truck_list)):length(used_trucks); 
    #     by=t_truck -> truck_sort_fn(t_truck, ntrucks, TR)
    # )

    # start = time()
    # # TODO too long?
    # for t_truck in candi_t_truck_list
    #     # push!(used_trucks, t_truck)
    #     # partialsort!(
    #     #     used_trucks, 
    #     #     length(used_trucks); 
    #     #     by=t_truck -> truck_sort_fn(t_truck, ntrucks, TR)
    #     # )
    #     insert_ind = findfirst(t_tr -> truck_sort_fn(t_truck, ntrucks, item_dispatch, TR, compatible_items) < truck_sort_fn(t_tr, ntrucks, item_dispatch, TR, compatible_items), used_trucks)
    #     insert!(
    #         used_trucks, 
    #         isnothing(insert_ind) ? lastindex(used_trucks)+1 : insert_ind, 
    #         t_truck
    #     )
    # end
    # println("Inserting: done in $(time() - start)s")

    # if !issorted(used_trucks, by=x -> truck_sort_fn(x, ntrucks, TR))
    #     display(candi_t_truck_list)
    #     # could fire sometimes if not sorted
    #     error("used_trucks is not sorted.")
    # end
end

function shufflesorted(seq; by=identity)
    res = []
    start = 1
    for (i, elem) in enumerate(seq)
        if i == 1
            continue
        end
        if by(elem) > by(seq[i-1])
            append!(res, shuffle(seq[start:i-1]))
            start = i
        end
    end
    append!(res, shuffle(seq[start:end]))
    return res
end

function solve_tsi_step!(i_items, item_dispatch, item_undispatch, used_trucks, solution, item_index; reshuffles=1)
    
    # notplaced_global = Pair{Int64, Item}[]
    # solve subproblems, retrieve items which couldn't be placed
    println()
    for (c, (t, truck)) in enumerate(used_trucks)
        items_dispatched = [i_items[ind][2] for ind in 1:length(item_dispatch) if item_dispatch[ind] == t]
        count_before = length(items_dispatched)
        if t in keys(solution)
            items_in_sol = [item for (i, stack) in solution[t] for item in get_items(stack)]
            if length(items_dispatched) == length(items_in_sol) && length(filter(x -> !(x in items_in_sol), items_dispatched)) == 0
                continue
            end
        end
        best = nothing
        for shuffle in 1:reshuffles
            sol, notplaced = BLtruck(items_dispatched, truck)
            if isnothing(best) || length(best[2]) > length(notplaced)
                best = sol, notplaced
            end
            if shuffle < reshuffles
                sort!(items_dispatched, by=item -> (
                    get_supplier_order(truck, get_supplier(item)), 
                    get_supplier_dock_order(truck, get_supplier(item), get_supplier_dock(item)), 
                    get_plant_dock_order(truck, get_plant_dock(item))))

                # shuffle!(items_dispatched)
                items_dispatched::Vector{Item} = shufflesorted(items_dispatched; by=item -> (
                    get_supplier_order(truck, get_supplier(item)), 
                    get_supplier_dock_order(truck, get_supplier(item), get_supplier_dock(item)), 
                    get_plant_dock_order(truck, get_plant_dock(item))))
            end
        end
        solution[t], notplaced = best
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
        id_cp = item -> string(get_id(item), "::", get_copy_number(item))
        notplaced_indices = [item_index[id_cp(item)][1] for item in notplaced]

        item_dispatch[notplaced_indices] .= -1
        for ind in notplaced_indices
            push!(item_undispatch[ind], t)
        end
        # clean up solution
        filter!((s_stack) -> 
            get_pos(s_stack[2]).x + get_dim(s_stack[2]).le <= get_dim(truck).le, solution[t]
        )
        
        # transfer notplaced to notplaced_global
        # append!(notplaced_global, map(item -> Pair(item_index[string(get_id(item), "::", get_copy_number(item))], item), notplaced))
        # remove notplaced items from item_dispatch of the truck
        # filter!(x -> !(string(get_id(x[2]), "::", get_copy_number(x[2])) in [string(get_id(y), "::", get_copy_number(y)) for y in notplaced]), item_dispatch[t])


        count_after = count(r -> r == t, item_dispatch) + length(notplaced)
        if count_before != count_after
            error("Items disappeared/duplicated\n$count_before != $count_after")
        end
        display_progress(c, length(used_trucks); name="Solving trucks")
    end
    


    
    # return notplaced_global

end

function select_multiple_trucks!(candi_list, item_dispatch, t_trucks, i_items, 
    items_permutation, ntrucks, nbuniqitems, truck_batch_size, TR, qte_fn, compatible_items)
    # add trucks as long as volume of candi_trucks < volume of i_items
    vol_trucks = 0.
    it_cursor = 1
    nbitems_left = count(t -> t == -1, item_dispatch)
    # vol_items = sum([get_volume(item) for (i, item) in i_items])
    vol_items = sum([qte_fn(i_items[ind][2]) for ind in items_permutation if item_dispatch[ind] == -1])

    println()
    while (!isnothing(truck_batch_size) || vol_trucks <= vol_items) && 
        (isnothing(truck_batch_size) || it_cursor <= truck_batch_size) && 
            it_cursor <= nbitems_left
        # TODO adding truck can take a very long time if too many to add
        # (double loop item/truck)

        candi_t, candi_truck = select_new_truck!(
            t_trucks, i_items[[ind for ind in items_permutation if item_dispatch[ind] == -1][it_cursor]], 
            ntrucks, nbuniqitems, TR, compatible_items, item_dispatch
        )
        push!(candi_list, Pair(candi_t, candi_truck))
        # item_dispatch[candi_t] = []
        if !isnothing(truck_batch_size)
            display_progress(it_cursor, truck_batch_size; name="Adding trucks")
        else
            display_progress(100*(vol_trucks/vol_items), 100; name="Adding trucks")
        end
        vol_trucks += qte_fn(candi_truck)
        it_cursor += 1
    end
    # sort!(t_trucks, by=x -> truck_sort_fn(x, ntrucks, TR)) # TODO is too long! (possibly)
    if !issorted(t_trucks, by=x -> truck_sort_fn(x, ntrucks, item_dispatch, TR, compatible_items))
        display([get_id(t_truck[2]) for t_truck in t_trucks])
        display(map(x -> truck_sort_fn(x, ntrucks, item_dispatch, TR, compatible_items), t_trucks))
        error("t_trucks is not sorted.")
    end
end

function solve_tsi(t_trucks, i_items, TR, assignment_fn; truck_batch_size=nothing, skipfirstpass=false, reshuffles=1)
    # i in i_items is line index of TR
    # t in t_trucks is column index of TR
    # extra trucks are multiples of original planned trucks in t_trucks

    nbuniqitems = length(Set(get_id(i_item[2]) for i_item in i_items))
    nbitems = length(i_items)
    item_index = Dict(string(get_id(i_items[ind][2]), "::", get_copy_number(i_items[ind][2])) => (ind, i_items[ind][1]) for ind in 1:nbitems)

    # println("Synchronized indices:")
    # println(
    #     setdiff(
    #         union(Set([i_item[1] for i_item in  i_items]), Set(collect(1:nbitems))), 
    #         intersect(union(Set([i_item[1] for i_item in  i_items]), Set(collect(1:nbitems))))
    #     )
    # )
    # readline()
    ntrucks = length(t_trucks)

    ###### DEBUG
    println("Computing...")
    mean_num_compatible_trucks = sum([sum(TR[:, ((i_item[1]-1) % nbuniqitems)+1]) for i_item in i_items])/nbitems

    println("Mean number of compatible trucks for each item: $(mean_num_compatible_trucks)")
    println("Number of trucks: $(ntrucks)")
    ###################
    # bind item index => list of planned truck indices
    compatible_trucks = Dict(i => [j for j in 1:ntrucks if TR[j, ((i-1) % nbuniqitems)+1]] for (i, item) in i_items)
    compatible_items = Dict(t => [j for j in 1:nbitems if TR[t, ((j-1) % nbuniqitems)+1]] for (t, truck) in t_trucks)
    # used_trucks = []
    used_trucks::Vector{Pair{Int64, Truck}} = copy(t_trucks)
    t_trucks = [Pair(t + ntrucks, truck) for (t, truck) in t_trucks]
    # item_dispatch = Dict{Any, Vector{Pair{Integer, Item}}}(t => [] for (t, truck) in used_trucks) # index is used trucks

    item_dispatch = Vector{Int64}(undef, nbitems) # -1 if not dispatch, else truck index
    # index of item_dispatch is 1:nbitems
    item_undispatch = Dict(ci => [] for ci in 1:nbitems) # list of truck indices for each item i_items index (not i)
    fill!(item_dispatch, -1)
    items_permutation = [ind for ind in 1:nbitems] # item_permutations stores indices of item_dispatch and i_items (same)

    solution = Dict()
    # sort trucks by increasing cost and by decreasing versatility
    sort!(t_trucks, by=t_truck -> truck_sort_fn(t_truck, ntrucks, item_dispatch, TR, compatible_items))
    sort!(used_trucks, by=t_truck -> truck_sort_fn(t_truck, ntrucks, item_dispatch, TR, compatible_items))

    # sort items by decreasing constraints
    # ((x[1]-1) % nbuniqitems)+1
    # sort!(i_items, by=i_item -> item_sort_fn(i_item, nbuniqitems, TR))
    # TODO sort items by similarity, so that similar objects are next to each other in the list
    sort!(items_permutation, by=ind -> item_sort_fn(Pair(i_items[ind][1], i_items[ind][2]), nbuniqitems, TR))

    candi_t, candi_truck = nothing, nothing
    println()
    display_progress(nbitems - count(t -> t != 1, item_dispatch) +  1, nbitems; name="Items placed")
    # while there are still items to place
    # while !isempty(i_items)
    while -1 in item_dispatch
        # DONE make select_new_truck add multiple trucks
        # DONE use intuition to find a lower bound on the number of trucks to add
        # DONE use intuition to find an upper bound on the number of trucks to add
        candi_list = Pair{Integer, Truck}[]
        if !skipfirstpass
            select_multiple_trucks!(
                candi_list, item_dispatch, t_trucks, i_items, items_permutation, ntrucks, 
                # nbuniqitems, truck_batch_size, TR, get_volume
                nbuniqitems, truck_batch_size, TR, get_area, compatible_items
                )
            # println("Truck selection added $(length(candi_list)) new trucks")
            # sort!(candi_list, by=t_truck -> truck_sort_fn(t_truck, ntrucks))
        else
            skipfirstpass=false 
        end
        # candi_list = isnothing(candi_t) ? [] : [Pair(candi_t, candi_truck)]
        assignment_fn((i_items, [ind for ind in 1:nbitems if item_dispatch[ind] == -1], 
            candi_list, used_trucks, TR, item_dispatch, ntrucks, 
            nbuniqitems, compatible_trucks, compatible_items, item_undispatch)
        )
        # assign_to_trucks!(
        #     i_items, [ind for ind in 1:nbitems if item_dispatch[ind] == -1], candi_list, used_trucks, TR, item_dispatch, ntrucks, 
        #     nbuniqitems, compatible_trucks; limit=item_batch_size
        # )
        solve_tsi_step!(i_items, item_dispatch, item_undispatch, used_trucks, solution, item_index; reshuffles=reshuffles)
        
        # Already done by design?
        if !issorted(used_trucks, by=x -> truck_sort_fn(x, ntrucks, item_dispatch, TR, compatible_items))

            error("used_trucks is not sorted.")
        end

        # for stuff in notplaced_global
        #     if stuff in i_items
        #         error("$stuff in i_items")
        #     end
        # end

        # append!(i_items, notplaced_global)

        # sort!(i_items, by= i_item -> item_sort_fn(i_item, nbuniqitems, TR)) # TODO insert in placed?
        
        clearnlines(2 + !skipfirstpass)
        display_progress(nbitems - count(t -> t == -1, item_dispatch) +  1, nbitems; name="Item placed")
    end
    
    return used_trucks, solution
end

function mk_benders_master(
    optimizer, coefcosttransportation::Real, 
    coefcostextratruck::Real, coefcostinventory::Real,
    t_trucks, _i_items, items_indices, _TR, IDL, TDA, IM, TMm, IV, TV, nbuniqitems, nbuniqtrucks,
    item_undispatch; silent=false, timeout=nothing
)
    i_items = _i_items[items_indices]
    nbitems = length(i_items)
    nbtrucks = length(t_trucks)

    TR = Matrix{Float64}(undef, nbtrucks, nbitems)
    for ci in items_indices
        i, item = _i_items[ci]
        for (ct, (t, truck)) in enumerate(t_trucks)
            TR[ct, ci] = t in item_undispatch[ci] ? 0. : _TR[((t-1) % nbuniqtrucks)+1, ((i-1) % nbuniqitems)+1]
        end
    end

    nbTR::Int64 = sum(TR)

    trucks = [t_truck[2] for t_truck in t_trucks]
    items = [i_item[2] for i_item in i_items]
    
# function TSIModel(optimizer, 
#     item_productcodes, truckdict, supplierdict, supplierdockdict, plantdict, 
#     plantdockdict, nbplannedtrucks, nbitems, nbsuppliers, nbsupplierdocks, nbplants, 
#     nbplantdocks, nbstacks, TE_P, TL_P, TW_P, TH_P, TKE_P, TGE_P, TDA_P, TU_P, TP_P, TK_P, 
#     TG_P, TR_P, IU, IP, IK, IPD, IS, IDL, IDE, stackabilitycodedict, nbtrucks, TE, 
#     TL, TW, TH, TKE, TGE, TDA, TU, TP, TK, TG, TR, TID, reverse_truckdict
#     )

    ## Create model
    # model = Model(
    #     optimizer_with_attributes(optimizer, "timeLimit" => timeout)
        
    #     )

    model = Model(
        optimizer
        )
        

    ## Add variables
    @info "Creating variables..."
    
    @info "Adding TI..."
    # @variable(model, TI[1:nbtrucks, 1:nbitems], lower_bound = 0, upper_bound = 1, Bin)
    @variable(model, TIvars[1:nbTR] >= 0)
    @constraint(model, TIvars .<= 1)
    
    TRindexes = findall(x -> x == 1, TR)
    TI = sparse([r[1] for r in TRindexes], [r[2] for r in TRindexes], TIvars)
    
    if size(TI)[1] < nbtrucks
        # println(size(TI))
        # println(size(zeros((nbtrucks - size(TI)[1]), size(TI)[2])))
        TI = sparse_vcat(TI, zeros((nbtrucks - size(TI)[1]), size(TI)[2]))
    end

    @info "Adding zeta..."
    @variable(model, zeta[1:nbtrucks] >= 0)
    @info "Adding zetab..."
    @variable(model, zetab[1:nbtrucks] >= 0)

    # display(TI)
    # readline()

    # TR_truck_index = Vector{Tuple{Integer, Integer}}(undef, nbtrucks) # index: truck index
    # fill!(TR_truck_index, (-1, -1))
    # # TR_truck_index[1] = (0, sum(TR[1, :]))
    # TR_item_index = Dict(ci => [] for ci in 1:nbitems) # maps item index to variables' indexes in TI
    # revTR_item_index = Dict() # maps variables' indexes in TI to item index
    # cnt = 1
    # for (ct, (t, truck)) in enumerate(t_trucks)
    #     if (-1, -1) in TR_truck_index
    #         before = ct-1 >= 1 ? TR_truck_index[ct-1][2] : 0
    #         TR_truck_index[ct] = before, before + sum(TR[ct, :])
    #     end
    #     for (ci, (i, item)) in enumerate(i_items)
    #         if TR[ct, ci] == 1
    #             push!(TR_item_index[ci], cnt)
    #             revTR_item_index[cnt] = ci
    #             cnt += 1
    #         end
    #     end
    # end

    # display(TR_truck_index)
    # readline()

    # println(findfirst(x -> sum(TR[:, x]) == 0, 1:nbitems))
    # readline()
    @info "Computing parameters..."

    Mzeta = 100000

    Alpha = diagm([get_id(trucks[t])[1] == 'P' ? coefcosttransportation : coefcostextratruck for t in 1:nbtrucks])
    # Alpha = diagm([get_id(trucks[t])[1] == 'P' ? coefcosttransportation : coefcostextratruck for t in 1:nbtrucks])

    ## Add constraints
    @info "Adding constraints..."
    @info "Adding czeta_TI_zetab..."
    # @constraint(model, czeta_TI_zetab, zeta .<= TI * ones(size(TI)[2]) - ones(size(zeta)[1]) + zetab)
    # println(size(zeta))
    # println(size(TI))
    # println(size(ones(size(TI)[2])))
    # println(size(ones(size(TI)[1])))
    # println(size(zetab))
    @constraint(model, czeta_TI_zetab, zeta .<= TI * ones(size(TI)[2]) - ones(size(TI)[1]) + zetab)

    @info "Adding czetab_TI..."
    # @constraint(model, czetab_TI, zetab .<= ones(size(zetab)[1]) - TI * ones(size(TI)[2]) * Mzeta)
    @constraint(model, czetab_TI_left, zetab .>= ones(size(TI)[1]) - TI * ones(size(TI)[2]) * Mzeta)
    @constraint(model, czetab_TI_right, zetab .<= ones(size(TI)[1]) + TI * ones(size(TI)[2]) * Mzeta)

    @info "Adding czetab_1_TI..."

    @constraint(model, czetab_1_TI, zetab .<= (ones(size(TI)[1]) - TI * ones(size(TI)[2])) * Mzeta)
    @info "Adding cTI_IM_TMm..."

    @constraint(model, cTI_IM_TMm, TI * IM .<= TMm)
    @info "Adding TI_IV_TV..."

    @constraint(model, TI_IV_TV, TI * IV .<= TV)
    @info "Adding TI_1_1..."

    @constraint(model, TI_1_1, transpose(TI) * ones(size(TI)[1]) == ones(size(TI)[2]))

    # @info "Adding cTI_TR..."
    # @constraint(model, cTI_TR, TI .<= TR)
    
    # DONE add global weight constraint
    # @info "Adding cTI_IM_TMm..."
    # @constraint(model, cTI_IM_TMm, TI * IM <= TMm)
    
    # @info "Adding cTI_1_1..."
    # display(size(transpose(TI)))
    # display(size(ones(size(transpose(TI))[2])))
    # display(size(transpose(TI) * ones(size(transpose(TI))[2])))
    # @constraint(model, cTI_1_1, transpose(TI) * ones(size(transpose(TI))[2]) == ones(size(transpose(TI))[1]))
    
    # zeta_cnt = 0
    # for (ct, (t, truck)) in enumerate(t_trucks)
    #     display_progress(ct, nbtrucks; name="Building Constraints (trucks)", margin=40)
    #     if sum(TR[ct, :]) == 0
    #         continue
    #     end
    #     zeta_cnt += 1
    #     TIct_it_indices = [i for i in 1:nbTR if TR_truck_index[ct][1] < i <= TR_truck_index[ct][2]]
    #     # println(TIct_it_indices)
    #     # @constraint(
    #     #     model, 
    #     #     zeta[zeta_cnt] <= transpose(TI[TIct_it_indices]) * ones(convert(Int64, sum(TR[ct, :]))) - 1 + zetab[zeta_cnt], 
    #     #     set_string_name=false # faster building in theory
    #     #     # base_name=string("czeta_TI_zetab", ct)
    #     # )

    #     # @constraint(
    #     #     model,  
    #     #     # zetab[ct] <= 1 - TI[TIct_it_indices] * ones(convert(Int64, sum(TR[ct, :]))) * Mzeta,
    #     #     zetab[zeta_cnt] >= 1 - sum(TI[TIct_it_indices]) * Mzeta,
    #     #     set_string_name=false # faster building in theory
    #     #     # base_name=string("czetab_TI_left", ct)
    #     # )

    #     # @constraint(
    #     #     model,  
    #     #     # zetab[ct] <= 1 - TI[TIct_it_indices] * ones(convert(Int64, sum(TR[ct, :]))) * Mzeta,
    #     #     zetab[zeta_cnt] <= 1 + sum(TI[TIct_it_indices]) * Mzeta,
    #     #     set_string_name=false # faster building in theory
    #     #     # base_name=string("czetab_TI_right", ct)
    #     # )
        
    #     # @constraint(
    #     #     model,  
    #     #     # zetab[ct] <= 1 - TI[TIct_it_indices] * ones(convert(Int64, sum(TR[ct, :]))) * Mzeta,
    #     #     zetab[zeta_cnt] <= (1 - sum(TI[TIct_it_indices])) * Mzeta,
    #     #     set_string_name=false # faster building in theory
    #     #     # base_name=string("czetab_1_TI", ct)
    #     # )
        

    #     # @constraint(
    #     #     model, 
    #     #     transpose(TI[TIct_it_indices]) * IM[[revTR_item_index[it] for it in TIct_it_indices]] <= TMm[ct],
    #     #     set_string_name=false # faster building in theory
    #     #     # base_name=string("cTI_IM_TMm", ct)    
    #     # )

    #     # println(size(transpose(TI[TIct_it_indices])))
    #     # println(size([get_volume(i_items[revTR_item_index[it]][2]) for it in TIct_it_indices]))
    #     # println(size(transpose(TI[TIct_it_indices]) * [get_volume(i_items[revTR_item_index[it]][2]) for it in TIct_it_indices]))
    #     # println(size(get_volume(t_trucks[ct][2])))

    #     # @constraint(
    #     #     model, 
    #     #     transpose(TI[TIct_it_indices]) * [get_volume(i_items[revTR_item_index[it]][2]) for it in TIct_it_indices] <= get_volume(t_trucks[ct][2]),
    #     #     set_string_name=false # faster building in theory
    #     #     # base_name=string("cTI_IV_TV", ct)    
    #     # )

    # end
    # for (ci, (i, item)) in enumerate(i_items)
    #     display_progress(ci, nbitems; name="Building Constraints (items)", margin=40)

    #     @constraint(
    #         model, 
    #         transpose(TI[TR_item_index[ci]]) * ones(length(TR_item_index[ci])) == 1,
    #         set_string_name=false # faster building in theory
    #         # base_name=string("cTI_1_1_", ci)
    #     )
    # end
        

    # Cost of used planned trucks
    @expression(model, obj1, sum(transpose([get_cost(trucks[t]) for t in 1:nbtrucks]) * Alpha * TI * ones(size(TI)[2])))
    # @expression(model, obj1, 
    #     sum(
    #         get_cost(trucks[ct]) * Alpha[numt, numt] * sum(TI[[i for i in 1:nbTR if TR_truck_index[ct][1] < i <= TR_truck_index[ct][2]]])
    #         for (numt, ct) in enumerate([t for t in 1:nbtrucks if 1 in TR[t, :]])
    #     )
    # )

    # Cost of used extra trucks
    @expression(model, obj2, 
        -sum(
            transpose(zeta) * Alpha * [get_cost(trucks[t]) for t in 1:nbtrucks]
        )
    
    )
    # inventory cost
    @expression(model, obj3, coefcostinventory * transpose([get_inventory_cost(items[i]) for i in 1:nbitems]) * (IDL - transpose(TI) * TDA))
    
    # display(findfirst(cj -> TR_truck_index[ct][1] < cj <= TR_truck_index[ct][2], TR_item_index[1]))
    # display(!isnothing(findfirst(cj -> TR_truck_index[ct][1] < cj <= TR_truck_index[ct][2], TR_item_index[1])))
    # println(size(TDA[[
    #     ct for ct in 1:nbtrucks if !isnothing(findfirst(cj -> TR_truck_index[ct][1] < cj <= TR_truck_index[ct][2], TR_item_index[1]))
    # ]]))
    # println(sum(TR[:, 1]))
    # println(size(transpose(TR[TR_item_index[1]])))

    # println(size(transpose([get_inventory_cost(items[i]) for i in 1:nbitems])))
    # println(size(IDL))
    # @expression(model, obj3, 
    #     coefcostinventory * 
    #     transpose([get_inventory_cost(items[i]) for i in 1:nbitems]) * 
    #     (
    #         IDL - [
    #             transpose(TR[TR_item_index[ci]]) * 
    #             # TDA[[
    #             #     ct for ct in 1:nbtrucks if !isnothing(findfirst(cj -> TR_truck_index[ct][1] < cj <= TR_truck_index[ct][2], TR_item_index[ci]))
    #             # ]]
    #             TDA[findall(ct -> TR[ct, ci] == 1, 1:nbtrucks)]
    #             for ci in 1:nbitems
    #         ]
    #     )
    # )
    # @objective(submodel, Min, 
    # sum(subpb[:costtransportation] * submodel[:zetaT]) + 
    # sum(subpb[:costextratruck] * submodel[:zetaE]) + 
    # subpb[:costinventory] * sum(subpb[:IDL] - submodel[:TI][t, :] * subpb[:TDA][t]) + 
    # kappa * sum(submodel[:TI] - TIbar))

    @objective(model, Min, obj1 + obj2 + obj3)
    if silent
        set_silent(model)
    end

    # return model, TR_truck_index, TR_item_index, revTR_item_index
    return model, nothing, nothing, nothing
end

function solve_benders_master(
    model; relax=false, timeout=nothing
)
    if !isnothing(timeout)
        set_time_limit_sec(model, timeout)
    end
    if relax
        undo = relax_integrality(model)
    end
    optimize!(model)
    # TIvalues = Matrix{Float64}(undef, size(model[:TI])...)
    TIvalues = Vector{Float64}(undef, size(model[:TI])...)
    # TIvalues .= round.(value.(model[:TI])) # TODO round is risky?
    TIvalues .= value.(model[:TI]) # TODO round is risky?

    return TIvalues
end

# assign_benders!(::Vector{Pair{Integer, Item}}, ::Vector{Int64}, ::Vector{Pair{Integer, Truck}}, ::Vector{Pair{Int64, Truck}}, ::BitMatrix, ::Vector{Int64}, ::Int64, ::Int64, ::Dict{Int64, Vector{Int64}}, ::Dict{Int64, Vector{Any}}, ::Type{Cbc.Optimizer}, ::Float64, ::Float64, ::Float64; relax=true, silent=false)
# assign_benders!(::Vector{Pair{Integer, Item}}, ::Vector{<:Integer}, ::Vector{Pair{Integer, Truck}}, ::Vector{Pair{Integer, Truck}}!!!, ::Any, ::Vector{<:Integer}, ::Integer, ::Integer, ::Any, ::Any, ::Any, ::Real, ::Real, ::Real; relax, silent)

function assign_benders!(i_items::Vector{Pair{Integer, Item}}, items_indices::Vector{<:Integer}, 
    candi_t_truck_list::Vector{Pair{Integer, Truck}}, used_trucks::Vector{Pair{T, Truck}}, 
    TR, item_dispatch::Vector{<:Integer}, nbuniqtrucks::Integer, 
    nbuniqitems::Integer, compatible_trucks, item_undispatch, optimizer, coefcosttransportation::Real, 
    coefcostextratruck::Real, coefcostinventory::Real; relax=false, silent=false, timeout=nothing) where T <: Integer

    #### DEBUG
    tmp_trucklist = [p for truck_list in [candi_t_truck_list, used_trucks] for p in truck_list]
    # append!(tmp_trucklist, tmp_trucklist)
    # append!(tmp_trucklist, tmp_trucklist)
    # append!(tmp_trucklist, tmp_trucklist)
    # ((i_items[ind][1]-1) % nbuniqitems)+1
    # items_indices = [
    #     ind for ind in items_indices if sum(TR[[p[1] for p in tmp_trucklist], ((i_items[ind][1]-1) % nbuniqitems)+1]) >= 1
    # ]
    ###################################################
    
    # TODO find upper bound for number of trucks

    TMm = [get_TMm(truck) for (t, truck) in tmp_trucklist]
    TV = [get_volume(truck) for (t, truck) in tmp_trucklist]

    IM = [get_weight(item) for (i, item) in i_items[items_indices]]
    IV = [get_volume(item) for (i, item) in i_items[items_indices]]

    TDA = [get_arrival_time(truck) for (t, truck) in tmp_trucklist]

    # converted to days
    IDL = [to_days(string(convert(Int64, get_time_window(item).latest))) for (i, item) in i_items[items_indices]]


    model, TR_truck_index, TR_item_index, revTR_item_index = mk_benders_master(
        optimizer, coefcosttransportation, coefcostextratruck, coefcostinventory,
        tmp_trucklist, i_items, items_indices, TR, 
        IDL, TDA, IM, TMm, IV, TV, nbuniqitems, nbuniqtrucks, item_undispatch; silent=silent, timeout=timeout
    )
    write_to_file(model, "model.lp")
    # print(model)
    TIvalues = solve_benders_master(model; relax=relax, timeout=timeout)

    println("TIVALUES")
    display(sum(TIvalues))
    println(TIvalues)
    readline()

    for (ci, ind) in enumerate(items_indices)
        i_item = i_items[ind]
        # t = findfirst(x -> x == 1, TIvalues[:, ci])
        # if !isnothing(t)
        #     item_dispatch[ind] = t
        # end

        cj = TR_item_index[ci][findfirst(x -> TIvalues[x] == 1, TR_item_index[ci])]
        ct = findfirst(cr -> TR_truck_index[cr][1] < cj <= TR_truck_index[cr][2], 1:length(tmp_trucklist))
        item_dispatch[ind] = tmp_trucklist[ct][1]
    end
    println("item_dispatch")
    println(item_dispatch)
    readline()

end