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
                # result["max_stackability"][i],
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
                        result["max_stackability"][i],
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
    # println("Max stackability of 0090014400_03072022000875 (before solving)")
    # println(get_max_stackability(i_items[findfirst(x -> get_id(x[2]) == "0090014400_03072022000875", i_items)][2]))
    # readline()

    # println("Max stackability of 0090014400_03072022000875's product (before solving)")
    # println(get_max_stackability(get_product(i_items[findfirst(x -> get_id(x[2]) == "0090014400_03072022000875", i_items)][2])))
    # readline()

    # println("Pbtic Product")
    # println(get_product(i_items[findfirst(x -> get_id(x[2]) == "0090014400_03072022000875", i_items)][2]))
    # readline()
    ###
    coeftransportationcost, coefcostextratruck = result["costtransportation"], result["coefcostextratruck"]
    ###
    # item_batch_size = convert(Int64, ceil(length(i_items)/10))
    item_batch_size = nothing
    skipfirstpass=true
    assign_fn = x -> assign_to_trucks!(x...; limit=item_batch_size)
    # truck_batch_size = convert(Int64, ceil(sqrt(length(i_items))))
    truck_batch_size = convert(Int64, ceil(length(i_items)/10))
    # truck_batch_size = nothing
    ###
    reshuffles = 2
    fusion_used_trucks, fusion_solution, used_trucks, solution = solve_tsi(t_trucks, i_items, result["TR_P"], assign_fn, coeftransportationcost, coefcostextratruck; 
        truck_batch_size=truck_batch_size, skipfirstpass=skipfirstpass, reshuffles=reshuffles)
    
    
    println("Solved $(length(i_items_origin)) items in $(time() - start)s")

    nb_items_after = sum([
        length(get_items(s_stack[2])) 
        for t in keys(solution) 
        for s_stack in solution[t]
    ])
    println("Items after solving: $(nb_items_after)")



    len_used_trucks_before = length(used_trucks)

    # remove empty trucks from solution
    filter!(t_stacks -> !isempty(t_stacks[2]), solution)
    print("Solution associated to non-existent truck: ")
    println(findfirst(t -> isnothing(findfirst(t_truck -> t_truck[1] == t, used_trucks)), collect(keys(solution))))

    # remove empty trucks from used trucks 
    filter!(t_truck -> t_truck[1] in keys(solution), used_trucks)

    len_used_trucks_after = length(used_trucks)

    println("Removed $(len_used_trucks_before - len_used_trucks_after) empty trucks from solution")

    println("Number of trucks in solution: $(len_used_trucks_after)")
    println("Number of trucks in fusion solution: $(length(fusion_used_trucks))")
        
    for (sol, ut) in [(solution, used_trucks), (fusion_solution, fusion_used_trucks)]

        filling_rates = [
            (isempty(sol[t]) ? 0. : get_loaded_volume([stack for (i, stack) in sol[t]])) / get_volume(Dict(ut)[t]) 
                for t in keys(sol)
        ]

        # println("Max stackability of 0090014400_03072022000875")
        # println(get_max_stackability(i_items[findfirst(x -> get_id(x[2]) == "0090014400_03072022000875", i_items)][2]))
        # readline()
        # for t in keys(solution)
        #     for (si, s) in solution[t]
        #         if si == 2 && get_id(used_trucks[findfirst(x -> x[1] == t, used_trucks)][2]) == "P202273201"
        #             display(min([get_max_stackability(item) for item in get_items(s)]...))
        #             display(get_minmax_stackability(s))
        #             readline()
        #         end
        #         if min([get_max_stackability(item) for item in get_items(s)]...) != get_minmax_stackability(s)
        #             readline()
        #         end
        #     end
        # end

        println("Filling rate stats over $(length(ut)) trucks")
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


        nb_stacks_truck = [
            isempty(sol[t]) ? 0 : length(sol[t]) 
                for t in keys(sol)
        ]

        println("Number of stacks stats over $(length(ut)) trucks")
        println("Mean number of stacks per truck")
        display(mean(nb_stacks_truck))
        println("Min number of stacks per truck")
        display(min(nb_stacks_truck...))
        println("25% quantile number of stacks per truck")
        display(quantile(nb_stacks_truck, 0.25))
        println("Median number of stacks per truck")
        display(median(nb_stacks_truck))
        println("75% quantile number of stacks per truck")
        display(quantile(nb_stacks_truck, 0.75))
        println("Max number of stacks per truck")
        display(max(nb_stacks_truck...))
    end
    used_truck_dict = Dict(fusion_used_trucks)

    # display(solution)
    output_directory = instancepath
    if !isdir(output_directory)
        mkdir(output_directory)
    end
    append=false
    for t in keys(fusion_solution)
        truck = used_truck_dict[t]
        # reassign ids because BLtruck creates new stacks id-less
        for (i, stack) in fusion_solution[t]
            set_id!(stack, string(get_id(truck), "_", i))
        end
        solution_to_csv(truck, fusion_solution[t], output_directory, append=append)
        if !append
            append = true
        end
    end
end