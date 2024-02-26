using Test
using Statistics

using JuMP
using MathOptInterface

using Cbc
using GLPK
using Clp

const MOI = MathOptInterface


include("../src/assignment.jl")
include("../src/instance_loader.jl")
include("../src/to_csv.jl")

function get_docks(supplier, supplierdockdict)
    return [split(k, "__")[2] for k in keys(supplierdockdict) if split(k, "__")[1] == supplier]
end

@testset "solve_tsi" begin
    instancepath = "./instances/BY2/"
    start = time()
    result = loadinstance(instancepath; onlyplanned=true)
    println("Loaded instance in $(time() - start)s")
    supplierdict = result["supplierdict"]
    supplierdockdict = result["supplierdockdict"]
    plantdockdict = result["plantdockdict"]
    plantdict = result["plantdict"]
    itemdict = result["itemdict"]
    # display(size(result["TR_P"]))
    nbplannedtrucks = result["nbplannedtrucks"]
    nbitems = result["nbitems"]
    t_trucks = []
    # display(supplierdict)
    # display(supplierdockdict)
    # display(get_docks(collect(keys(supplierdict))[1], supplierdockdict))
    # display(result["TE_P"])
    println()
    for t in 1:nbplannedtrucks
        display_progress(t, nbplannedtrucks; name="Instantiate trucks")
        # display([result["reverse_truckdict"][t],
            # Dim(result["TL_P"][t], result["TW_P"][t]),
            # result["TH_P"][t],
            # result["TEM_P"][t],
            # result["max_stack_weights_P"][result["reverse_truckdict"][t]],
            # result["TMm_P"][t],
            # Dict(
            #     parse(Int64, s) => Int(result["TE_P"][t, supplierdict[s]]) 
            #         for s in keys(supplierdict) if result["TU_P"][t, supplierdict[s]]
            # ),
            # Dict(parse(Int64, s) => Dict(sd => Int(result["TKE_P"][t, supplierdockdict[string(s, "__", sd)]]) for sd in get_docks(s, supplierdockdict)) 
            #     for s in keys(supplierdict) if result["TU_P"][t, supplierdict[s]]),
            # Dict(pd => result["TGE_P"][t, plantdockdict[string(p, "__", pd)]] for p in keys(plantdict) for pd in get_docks(p, plantdockdict)
            #      if result["TP_P"][t, plantdict[p]]),
            # result["Cost_P"][t],
            # result["CM_P"][t],
            # result["CJfm_P"][t],
            # result["CJfc_P"][t],
            # result["CJfh_P"][t],
            # result["EM_P"][t],
            # result["EJhr_P"][t],
            # result["EJcr_P"][t],
            # result["EJeh_P"][t],
            # result["EMmr_P"][t],
            # result["EMmm_P"][t]]
        # )
        # readline()
        push!(
            t_trucks, 
            Pair(
                t, 
                Truck(
                    result["reverse_truckdict"][t],
                    Dim(result["TL_P"][t], result["TW_P"][t]),
                    result["TH_P"][t],
                    result["TEM_P"][t],
                    result["max_stack_weights_P"][result["reverse_truckdict"][t]],
                    result["TMm_P"][t],
                    Dict(
                        parse(Int64, s) => Int(result["TE_P"][t, supplierdict[s]]) 
                            for s in keys(supplierdict) if result["TU_P"][t, supplierdict[s]]
                    ),
                    Dict(parse(Int64, s) => Dict(sd => Int(result["TKE_P"][t, supplierdockdict[string(s, "__", sd)]]) for sd in get_docks(s, supplierdockdict)) 
                        for s in keys(supplierdict) if result["TU_P"][t, supplierdict[s]]),
                    Dict(pd => result["TGE_P"][t, plantdockdict[string(p, "__", pd)]] for p in keys(plantdict) for pd in get_docks(p, plantdockdict)
                         if result["TP_P"][t, plantdict[p]]),
                    result["TDA_P"][t],
                    result["Cost_P"][t],
                    result["CM_P"][t],
                    result["CJfm_P"][t],
                    result["CJfc_P"][t],
                    result["CJfh_P"][t],
                    result["EM_P"][t],
                    result["EJhr_P"][t],
                    result["EJcr_P"][t],
                    result["EJeh_P"][t],
                    result["EMmr_P"][t],
                    result["EMmm_P"][t]
                )
            )
        )
    end
    # DONE set max_stack_weights
    i_items = Pair{Integer, Item}[]
    # display(t_trucks)
    products = Dict()
    println()
    for i in 1:nbitems
        display_progress(i, nbitems; name="Instantiate items")
        # get true nb items
        nbcopies = result["item_copies"][i]
        if !(result["item_productcodes"][i] in keys(products))
            products[result["item_productcodes"][i]] = Product(
                result["item_productcodes"][i],
                result["max_stackability"][i],
                100000 # I am just realizing that this is the same parameter as max_stack_weight in trucks
            )
        end

        for j in 1:nbcopies
            plantdockdict_keys = collect(keys(plantdockdict))
            plant_dock_name = plantdockdict_keys[
                findfirst(
                    y -> plantdockdict[y] == findfirst(
                        x -> x == 1, result["IPD"][i, :]
                    ), 
                    plantdockdict_keys
                )
            ]
            supplierdockdict_keys = collect(keys(supplierdockdict))
            supplier_dock_name = supplierdockdict_keys[
                findfirst(
                    y -> supplierdockdict[y] == findfirst(
                        x -> x == 1, result["IK"][i, :]
                    ), 
                    supplierdockdict_keys
                )
            ]

            push!(i_items,
                Pair(
                    i + nbitems * (j - 1),
                    Item(
                        itemdict[i],
                        result["pkg_code"][i],
                        j,
                        (earliest=result["IDE"][i], latest=result["IDL"][i]),
                        Dim(result["IL"][i], result["IW"][i]),
                        result["IH"][i],
                        result["IM"][i],
                        result["stackability_code"][i],
                        ismissing(result["_IO"][i]) ? :none : result["_IO"][i] == 1 ? :widthwise : :lengthwise,
                        split(plant_dock_name, "__")[1],
                        split(plant_dock_name, "__")[2],
                        parse(Int64, split(supplier_dock_name, "__")[1]),
                        split(supplier_dock_name, "__")[2],
                        result["inventory_cost"][i],
                        result["InH"][i],
                        products[result["item_productcodes"][i]]
                    )
                )
            )
        end
    end

    i_items_origin = copy(i_items)
    start = time()
    # display(i_items)

    ###
    # # optimizer = GLPK.Optimizer # Cbc or GLPK
    # optimizer = Clp.Optimizer # Cbc or GLPK
    # skipfirstpass=false
    # timeout = 60*2 #timelimit does nothing with Cbc, but works with Clp
    # truck_batch_size=convert(Int64, ceil(length(i_items)/10)) # never finishes, not feasible?
    # # truck_batch_size=convert(Int64, ceil(length(i_items)/2)) # too many trucks, takes too long
    # assign_fn = x -> assign_benders!(
    #     x..., 
    #     optimizer, 
    #     result["costtransportation"], 
    #     result["coefcostextratruck"],
    #     result["coefcostinventory"]; 
    #     relax=true, silent=false, timeout=timeout #relax = true is necessary for reasonable solving times
    # )
    ###

    ###
    # item_batch_size = convert(Int64, ceil(length(i_items)/10))
    item_batch_size = nothing
    skipfirstpass=true
    assign_fn = x -> assign_to_trucks!(x...; limit=item_batch_size)
    # truck_batch_size = convert(Int64, ceil(sqrt(length(i_items))))
    # truck_batch_size = convert(Int64, ceil(length(i_items)/10))
    truck_batch_size = nothing
    ###
    reshuffles = 2

    used_trucks, solution = solve_tsi(t_trucks, i_items, result["TR_P"], assign_fn; 
        truck_batch_size=truck_batch_size, skipfirstpass=skipfirstpass, reshuffles=reshuffles)
    
    
    println("Solved $(length(i_items_origin)) items in $(time() - start)s")

    len_used_trucks_before = length(used_trucks)

    # remove empty trucks from solution
    filter!(t_stacks -> !isempty(t_stacks[2]), solution)

    # remove empty trucks from used trucks 
    filter!(t_truck -> t_truck[1] in keys(solution), used_trucks)

    len_used_trucks_after = length(used_trucks)

    println("Removed $(len_used_trucks_before - len_used_trucks_after) empty trucks")

    filling_rates = [
        (isempty(solution[t]) ? 0. : get_loaded_volume([stack for (i, stack) in solution[t]])) / get_volume(Dict(used_trucks)[t]) 
            for t in keys(solution)
    ]

    println("Filling rate stats over $(length(used_trucks)) trucks")
    println("Mean filling rate")
    display(mean(filling_rates))
    println("Min filling rate")
    display(min(filling_rates...))
    println("25% quantile filling rate")
    display(quantile(filling_rates, 0.25))
    println("Median filling rate")
    display(median(filling_rates))
    println("75% quantile filling rate")
    display(quantile(filling_rates, 0.75))
    println("Max filling rate")
    display(max(filling_rates...))


    nb_items_truck = [
        isempty(solution[t]) ? 0 : length(solution[t]) 
            for t in keys(solution)
    ]

    println("Number of stacks stats over $(length(used_trucks)) trucks")
    println("Mean number of stacks per truck")
    display(mean(nb_items_truck))
    println("Min number of stacks per truck")
    display(min(nb_items_truck...))
    println("25% quantile number of stacks per truck")
    display(quantile(nb_items_truck, 0.25))
    println("Median number of stacks per truck")
    display(median(nb_items_truck))
    println("75% quantile number of stacks per truck")
    display(quantile(nb_items_truck, 0.75))
    println("Max number of stacks per truck")
    display(max(nb_items_truck...))
    
    used_truck_dict = Dict(used_trucks)

    # display(solution)
    output_directory = instancepath
    if !isdir(output_directory)
        mkdir(output_directory)
    end
    append=false
    for t in keys(solution)
        truck = used_truck_dict[t]
        # reassign ids because BLtruck creates new stacks id-less
        for (i, stack) in solution[t]
            set_id!(stack, string(get_id(truck), "_", i))
        end
        solution_to_csv(truck, solution[t], output_directory, append=append)
        if !append
            append = true
        end
    end
end