using JuMP
using LinearAlgebra

using Dates

include("dim.jl")
include("truck.jl")
include("item.jl")
include("progress.jl")
include("placement_algorithms.jl")
include("placement_visualizer.jl")

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

truck_sort_fn(t_truck, ntrucks) = (t_truck[1] <= ntrucks ? 0. : get_cost(t_truck[2]), ((.*).(versatility(t_truck[2]), -1))...)

item_sort_fn(i_item, nbuniqitems, TR) = (.*).(constraints(i_item[2], ((i_item[1]-1) % nbuniqitems)+1, TR), -1)

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

function select_new_truck!(t_trucks, i_item, ntrucks, nbuniqitems, TR)
    ## New truck phase
    # take cheapest truck that can take first item (with most constraints)
    candi_t, candi_truck = t_trucks[findfirst(
        t_truck -> valid_truck(t_truck, i_item, ntrucks, nbuniqitems, 0., 0., TR),
        t_trucks
    )]
    # remove the truck candidate from candidate trucks
    filter!(t_truck -> t_truck[1] != candi_t, t_trucks)
    
    # candi_truck = copy(candi_truck)
    candi_truck = set_id(candi_truck, new_truck_id(candi_t, candi_truck, ntrucks)) # change P to Q if extra + _ + copy number
    
    # Add extra truck to candidate trucks
    push!(t_trucks, Pair(candi_t + ntrucks, candi_truck)) # TODO keep it sorted
    sort!(t_trucks, by=x -> truck_sort_fn(x, ntrucks))
    if !issorted(t_trucks, by=x -> truck_sort_fn(x, ntrucks))
        display([get_id(t_truck[2]) for t_truck in t_trucks])
        display(map(x -> truck_sort_fn(x, ntrucks), t_trucks))
        error("t_trucks is not sorted.")
    end
    

    return candi_t, candi_truck

end

function assign_to_trucks!(
    i_items, items_indices, candi_t_truck_list, used_trucks, TR, item_dispatch, ntrucks, 
    nbuniqitems, compatible_trucks; limit=nothing
)

    tmp_trucklist = [p for truck_list in [candi_t_truck_list, used_trucks] for p in truck_list]

    items_volume = Dict(t => 0. for (t, truck) in tmp_trucklist)
    items_weight = Dict(t => 0. for (t, truck) in tmp_trucklist)

    # TODO function too long to process
    # TODO Idea: keep in memory a list of truck indices for each "group of items" for which those trucks work
    # Then, only probe those trucks for the items instead of trying everything
    ## Assigning to trucks phase
    # then for each item remaining to place
    for (c, ind) in enumerate(items_indices)
        i, item = i_items[ind]
        if !isnothing(limit) && c > limit
            break
        end
        for (t, truck) in tmp_trucklist
            if !(((t-1) % ntrucks)+1 in compatible_trucks[((i-1) % nbuniqitems)+1])
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
    append!(used_trucks, candi_t_truck_list) # TODO keep used_trucks sorted
    sort!(used_trucks, by=t_truck -> truck_sort_fn(t_truck, ntrucks))
    # if !issorted(used_trucks, by=x -> truck_sort_fn(x, ntrucks))
    #     display(candi_t_truck_list)
    #     # could fire sometimes if not sorted
    #     error("used_trucks is not sorted.")
    # end
end

function solve_tsi_step!(i_items, item_dispatch, used_trucks, solution, item_index)
    
    # notplaced_global = Pair{Int64, Item}[]
    # solve subproblems, retrieve items which couldn't be placed
    for (c, (t, truck)) in enumerate(used_trucks)
        items_dispatched = [i_items[ind][2] for ind in 1:length(item_dispatch) if item_dispatch[ind] == t]
        count_before = length(items_dispatched)
        solution[t], notplaced = BLtruck(items_dispatched, truck)
        # remove items of stacks which poke out of truck
        append!(
            notplaced, 
            [
                item 
                for (s, stack) in solution[t] 
                for item in get_items(stack)
                if get_pos(stack).x + get_dim(stack).le > get_dim(truck).le # when commented no more pb?
            ]
        )
        item_dispatch[[item_index[string(get_id(item), "::", get_copy_number(item))][1] for item in notplaced]] .= -1
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

function solve_tsi(t_trucks, i_items, TR, assignment_fn; truck_batch_size=nothing)
    # i in i_items is line index of TR
    # t in t_trucks is column index of TR
    # extra trucks are multiples of original planned trucks in t_trucks

    nbuniqitems = length(Set(get_id(i_item[2]) for i_item in i_items))
    nbitems = length(i_items)
    item_index = Dict(string(get_id(i_items[ind][2]), "::", get_copy_number(i_items[ind][2])) => (ind, i_items[ind][1]) for ind in 1:nbitems)

    ntrucks = length(t_trucks)

    ###### DEBUG
    println("Computing...")
    mean_num_compatible_trucks = sum([sum(TR[:, ((i_item[1]-1) % nbuniqitems)+1]) for i_item in i_items])/nbitems

    println("Mean number of compatible trucks for each item: $(mean_num_compatible_trucks)")
    println("Number of trucks: $(ntrucks)")
    ###################
    # bind item index => list of planned truck indices
    compatible_trucks = Dict(i => [j for j in 1:ntrucks if TR[j, ((i-1) % nbuniqitems)+1]] for (i, item) in i_items)

    # used_trucks = []
    used_trucks::Vector{Pair{Integer, Truck}} = copy(t_trucks)
    t_trucks = [Pair(t + ntrucks, truck) for (t, truck) in t_trucks]
    # item_dispatch = Dict{Any, Vector{Pair{Integer, Item}}}(t => [] for (t, truck) in used_trucks) # index is used trucks

    item_dispatch = Vector{Int64}(undef, nbitems) # -1 if not dispatch, else truck index
    fill!(item_dispatch, -1)
    items_permutation = [ind for ind in 1:nbitems] # item_permutations stores indices of item_dispatch and i_items (same)

    solution = Dict()
    # sort trucks by increasing cost and by decreasing versatility
    sort!(t_trucks, by=t_truck -> truck_sort_fn(t_truck, ntrucks))
    sort!(used_trucks, by=t_truck -> truck_sort_fn(t_truck, ntrucks))

    # sort items by decreasing constraints
    # ((x[1]-1) % nbuniqitems)+1
    # sort!(i_items, by=i_item -> item_sort_fn(i_item, nbuniqitems, TR))
    sort!(items_permutation, by=ind -> item_sort_fn(Pair(i_items[ind][1], i_items[ind][2]), nbuniqitems, TR))

    firstpass = true
    candi_t, candi_truck = nothing, nothing
    display_progress(nbitems - count(t -> t != 1, item_dispatch) +  1, nbitems; name="Items placed")
    # while there are still items to place
    # while !isempty(i_items)
    while -1 in item_dispatch
            # println("Items left: $(length(i_items))")
        
        # println([iit[1] for iit in i_items])
        # display([length(item_dispatch[k]) for k in keys(item_dispatch)])
        # # display([iit[1] for k in keys(item_dispatch) for iit in item_dispatch[k]])
        # println("Sum of items = $(length(i_items) + sum([length(item_dispatch[k]) for k in keys(item_dispatch)]))")
        # readline()
        # error("Number of items increases")
        # TODO make select_new_truck add multiple trucks
        # TODO use intuition to find a lower bound on the number of trucks to add
        candi_list = Pair{Integer, Truck}[]
        if !firstpass
            # add trucks as long as volume of candi_trucks < volume of i_items
            vol_trucks = 0.
            it_cursor = 1
            nbitems_left = count(t -> t == -1, item_dispatch)
            # vol_items = sum([get_volume(item) for (i, item) in i_items])
            vol_items = sum([get_volume(i_items[ind][2]) for ind in items_permutation if item_dispatch[ind] == -1])
            if isnothing(truck_batch_size)
                println()
            end
            while (!isnothing(truck_batch_size) || vol_trucks < vol_items) && 
                (isnothing(truck_batch_size) || it_cursor <= truck_batch_size) && 
                    it_cursor <= nbitems_left
                # TODO adding truck can take a very long time if too many to add
                # (double loop item/truck)

                candi_t, candi_truck = select_new_truck!(
                    t_trucks, i_items[[ind for ind in items_permutation if item_dispatch[ind] == -1][it_cursor]], ntrucks, nbuniqitems, TR
                )
                push!(candi_list, Pair(candi_t, candi_truck))
                # item_dispatch[candi_t] = []
                if !isnothing(truck_batch_size)
                    display_progress(it_cursor, truck_batch_size; name="Adding trucks")
                else
                    display_progress(100*(vol_trucks/vol_items), 100; name="Adding trucks")
                end
                vol_trucks += get_volume(candi_truck)
                it_cursor += 1
            end
            # sort!(candi_list, by=t_truck -> truck_sort_fn(t_truck, ntrucks))
        else
            firstpass=false 
        end
        # candi_list = isnothing(candi_t) ? [] : [Pair(candi_t, candi_truck)]
        assignment_fn((i_items, [ind for ind in 1:nbitems if item_dispatch[ind] == -1], 
            candi_list, used_trucks, TR, item_dispatch, ntrucks, 
            nbuniqitems, compatible_trucks)
        )
        # assign_to_trucks!(
        #     i_items, [ind for ind in 1:nbitems if item_dispatch[ind] == -1], candi_list, used_trucks, TR, item_dispatch, ntrucks, 
        #     nbuniqitems, compatible_trucks; limit=item_batch_size
        # )
        solve_tsi_step!(i_items, item_dispatch, used_trucks, solution, item_index)
        
        # Already done by design?
        if !issorted(used_trucks, by=x -> truck_sort_fn(x, ntrucks))

            error("used_trucks is not sorted.")
        end

        # for stuff in notplaced_global
        #     if stuff in i_items
        #         error("$stuff in i_items")
        #     end
        # end

        # append!(i_items, notplaced_global)

        # sort!(i_items, by= i_item -> item_sort_fn(i_item, nbuniqitems, TR)) # TODO insert in placed?
        
        clearnlines(2 + !firstpass)
        display_progress(nbitems - count(t -> t == -1, item_dispatch) +  1, nbitems; name="Item placed")
    end
    
    return used_trucks, solution
end

function mk_benders_master(
    optimizer, coefcosttransportation::Real, 
    coefcostextratruck::Real, coefcostinventory::Real,
    t_trucks, i_items, _TR, IDL, TDA, IM, TMm, nbuniqitems, nbuniqtrucks; silent=false
)
    nbitems = length(i_items)
    nbtrucks = length(t_trucks)

    TR = Matrix{Float64}(undef, nbtrucks, nbitems)
    for (ci, (i, item)) in enumerate(i_items)
        for (ct, (t, truck)) in enumerate(t_trucks)
            TR[ct, ci] = _TR[((t-1) % nbuniqtrucks)+1, ((i-1) % nbuniqitems)+1]
        end
    end

    nbTR::Integer = sum(TR)

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
    model = Model(optimizer)

    ## Add variables
    @info "Creating variables..."
    @info "Adding zeta..."
    @variable(model, zeta[1:count(x -> 1 in TR[x, :], 1:nbtrucks)] >= 0)
    @info "Adding zetab..."
    @variable(model, zetab[1:count(x -> 1 in TR[x, :], 1:nbtrucks)] >= 0)

    @info "Adding TI..."
    # @variable(model, TI[1:nbtrucks, 1:nbitems], lower_bound = 0, upper_bound = 1, Bin)
    @variable(model, TI[1:nbTR], lower_bound = 0, upper_bound = 1, Bin)

    TR_truck_index = Vector{Tuple{Integer, Integer}}(undef, nbtrucks) # index: truck index
    fill!(TR_truck_index, (-1, -1))
    # TR_truck_index[1] = (0, sum(TR[1, :]))
    TR_item_index = Dict(ci => [] for ci in 1:nbitems) # maps item index to variables' indexes in TI
    revTR_item_index = Dict() # maps variables' indexes in TI to item index
    cnt = 1
    for (ct, (t, truck)) in enumerate(t_trucks)
        if (-1, -1) in TR_truck_index
            before = ct-1 >= 1 ? TR_truck_index[ct-1][2] : 0
            TR_truck_index[ct] = before, before + sum(TR[ct, :])
        end
        for (ci, (i, item)) in enumerate(i_items)
            if TR[ct, ci] == 1
                push!(TR_item_index[ci], cnt)
                revTR_item_index[cnt] = ci
                cnt += 1
            end
        end
    end

    # display(TR_truck_index)
    # readline()

    # println(findfirst(x -> sum(TR[:, x]) == 0, 1:nbitems))
    # readline()
    @info "Computing parameters..."

    Mzeta = 100000

    Alpha = diagm([get_id(trucks[t])[1] == 'P' ? coefcosttransportation : coefcostextratruck for t in 1:nbtrucks if 1 in TR[t, :]])
    # Alpha = diagm([get_id(trucks[t])[1] == 'P' ? coefcosttransportation : coefcostextratruck for t in 1:nbtrucks])

    ## Add constraints
    # @info "Adding constraints..."
    # @info "Adding czeta_TI_zetab..."
    # @constraint(model, czeta_TI_zetab, zeta .<= TI * ones(size(TI)[2]) - ones(size(zeta)[1]) + zetab)
    zeta_cnt = 0
    for (ct, (t, truck)) in enumerate(t_trucks)
        display_progress(ct, nbtrucks; name="Building Constraints")
        if sum(TR[ct, :]) == 0
            continue
        end
        zeta_cnt += 1
        TIct_it_indices = [i for i in 1:nbTR if TR_truck_index[ct][1] < i <= TR_truck_index[ct][2]]
        # println(TIct_it_indices)

        @constraint(
            model, 
            zeta[zeta_cnt] <= transpose(TI[TIct_it_indices]) * ones(convert(Int64, sum(TR[ct, :]))) - 1 + zetab[zeta_cnt], 
            base_name=string("czeta_TI_zetab", ct)
        )
        
        # @info "Adding czetab_TI..."
        # @constraint(model, czetab_TI, zetab .<= ones(size(zetab)[1]) - TI * ones(size(TI)[2]) * Mzeta)
        @constraint(
            model,  
            # zetab[ct] <= 1 - TI[TIct_it_indices] * ones(convert(Int64, sum(TR[ct, :]))) * Mzeta,
            zetab[zeta_cnt] >= 1 - sum(TI[TIct_it_indices]) * Mzeta,
            base_name=string("czetab_TI_left", ct)
        )

        @constraint(
            model,  
            # zetab[ct] <= 1 - TI[TIct_it_indices] * ones(convert(Int64, sum(TR[ct, :]))) * Mzeta,
            zetab[zeta_cnt] <= 1 + sum(TI[TIct_it_indices]) * Mzeta,
            base_name=string("czetab_TI_right", ct)
        )
        
        @constraint(
            model,  
            # zetab[ct] <= 1 - TI[TIct_it_indices] * ones(convert(Int64, sum(TR[ct, :]))) * Mzeta,
            zetab[zeta_cnt] <= (1 - sum(TI[TIct_it_indices])) * Mzeta,
            base_name=string("czetab_1_TI", ct)
        )
        

        # @info "Adding cTI_TR..."
        # @constraint(model, cTI_TR, TI .<= TR)
        
        # TODO add global weight constraint
        # @info "Adding cTI_IM_TMm..."
        # @constraint(model, cTI_IM_TMm, TI * IM <= TMm)
        @constraint(
            model, 
            transpose(TI[TIct_it_indices]) * IM[[revTR_item_index[it] for it in TIct_it_indices]] <= TMm[ct],
            base_name=string("cTI_IM_TMm", ct)    
        )

        # println(size(transpose(TI[TIct_it_indices])))
        # println(size([get_volume(i_items[revTR_item_index[it]][2]) for it in TIct_it_indices]))
        # println(size(transpose(TI[TIct_it_indices]) * [get_volume(i_items[revTR_item_index[it]][2]) for it in TIct_it_indices]))
        # println(size(get_volume(t_trucks[ct][2])))

        @constraint(
            model, 
            transpose(TI[TIct_it_indices]) * [get_volume(i_items[revTR_item_index[it]][2]) for it in TIct_it_indices] <= get_volume(t_trucks[ct][2]),
            base_name=string("cTI_IV_TV", ct)    
        )

    end
    # @info "Adding cTI_1_1..."
    # display(size(transpose(TI)))
    # display(size(ones(size(transpose(TI))[2])))
    # display(size(transpose(TI) * ones(size(transpose(TI))[2])))
    for (ci, (i, item)) in enumerate(i_items)
        # @constraint(model, cTI_1_1, transpose(TI) * ones(size(transpose(TI))[2]) == ones(size(transpose(TI))[1]))
        @constraint(
            model, 
            transpose(TI[TR_item_index[ci]]) * ones(length(TR_item_index[ci])) == 1,
            base_name=string("cTI_1_1_", ci)
        )
    end
        

    # Cost of used planned trucks
    # @expression(model, obj1, sum(transpose([get_cost(trucks[t]) for t in 1:nbtrucks]) * Alpha * TI * ones(size(TI)[2])))
    @expression(model, obj1, 
        sum(
            get_cost(trucks[ct]) * Alpha[numt, numt] * sum(TI[[i for i in 1:nbTR if TR_truck_index[ct][1] < i <= TR_truck_index[ct][2]]])
            for (numt, ct) in enumerate([t for t in 1:nbtrucks if 1 in TR[t, :]])
        )
    )

    # Cost of used extra trucks
    @expression(model, obj2, 
        -sum(
            transpose(zeta) * Alpha * [get_cost(trucks[t]) for t in 1:nbtrucks if 1 in TR[t, :]]
        )
    
    )
    # inventory cost
    # @expression(model, obj3, coefcostinventory * transpose([get_inventory_cost(items[i]) for i in 1:nbitems]) * (IDL - transpose(TI) * TDA))
    
    # display(findfirst(cj -> TR_truck_index[ct][1] < cj <= TR_truck_index[ct][2], TR_item_index[1]))
    # display(!isnothing(findfirst(cj -> TR_truck_index[ct][1] < cj <= TR_truck_index[ct][2], TR_item_index[1])))
    # println(size(TDA[[
    #     ct for ct in 1:nbtrucks if !isnothing(findfirst(cj -> TR_truck_index[ct][1] < cj <= TR_truck_index[ct][2], TR_item_index[1]))
    # ]]))
    # println(sum(TR[:, 1]))
    # println(size(transpose(TR[TR_item_index[1]])))

    # println(size(transpose([get_inventory_cost(items[i]) for i in 1:nbitems])))
    # println(size(IDL))
    @expression(model, obj3, 
        coefcostinventory * 
        transpose([get_inventory_cost(items[i]) for i in 1:nbitems]) * 
        (
            IDL - [
                transpose(TR[TR_item_index[ci]]) * 
                # TDA[[
                #     ct for ct in 1:nbtrucks if !isnothing(findfirst(cj -> TR_truck_index[ct][1] < cj <= TR_truck_index[ct][2], TR_item_index[ci]))
                # ]]
                TDA[findall(ct -> TR[ct, ci] == 1, 1:nbtrucks)]
                for ci in 1:nbitems
            ]
        )
    )
    # @objective(submodel, Min, 
    # sum(subpb[:costtransportation] * submodel[:zetaT]) + 
    # sum(subpb[:costextratruck] * submodel[:zetaE]) + 
    # subpb[:costinventory] * sum(subpb[:IDL] - submodel[:TI][t, :] * subpb[:TDA][t]) + 
    # kappa * sum(submodel[:TI] - TIbar))

    @objective(model, Min, obj1 + obj2 + obj3)
    if silent
        set_silent(model)
    end

    return model, TR_truck_index, TR_item_index, revTR_item_index
end

function solve_benders_master(
    model; relax=false
)
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



function assign_benders!(i_items::Vector{Pair{Integer, Item}}, items_indices::Vector{<:Integer}, 
    candi_t_truck_list::Vector{Pair{Integer, Truck}}, used_trucks::Vector{Pair{Integer, Truck}}, 
    TR, item_dispatch::Vector{<:Integer}, nbuniqtrucks::Integer, 
    nbuniqitems::Integer, compatible_trucks, optimizer, coefcosttransportation::Real, 
    coefcostextratruck::Real, coefcostinventory::Real; relax=false, silent=false)

    #### DEBUG
    tmp_trucklist = [p for truck_list in [candi_t_truck_list, used_trucks] for p in truck_list]
    append!(tmp_trucklist, tmp_trucklist)
    append!(tmp_trucklist, tmp_trucklist)
    append!(tmp_trucklist, tmp_trucklist)
    # ((i_items[ind][1]-1) % nbuniqitems)+1
    items_indices = [
        ind for ind in items_indices if sum(TR[[p[1] for p in tmp_trucklist], ((i_items[ind][1]-1) % nbuniqitems)+1]) >= 1
    ][end-100:end]
    ###################################################
    
    TMm = [get_TMm(truck) for (t, truck) in tmp_trucklist]

    IM = [get_weight(item) for (i, item) in i_items[items_indices]]

    TDA = [get_arrival_time(truck) for (t, truck) in tmp_trucklist]

    # converted to days
    IDL = [to_days(string(convert(Int64, get_time_window(item).latest))) for (i, item) in i_items[items_indices]]


    model, TR_truck_index, TR_item_index, revTR_item_index = mk_benders_master(
        optimizer, coefcosttransportation, coefcostextratruck, coefcostinventory,
        tmp_trucklist, i_items[items_indices], TR, 
        IDL, TDA, IM, TMm, nbuniqitems, nbuniqtrucks; silent=silent
    )
    write_to_file(model, "model.lp")
    # print(model)
    # TIvalues = solve_benders_master(model; relax=relax)
    TIvalues = solve_benders_master(model; relax=true) #TODO

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