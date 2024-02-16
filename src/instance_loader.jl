# module TSIInstanceLoader

# using LinearAlgebra
# using JuMP
# using Cbc
using CSV

# using FilePaths
# using Logging
# using Profile

# using ArgParse

function to_days(date::String)

    if '.' in date
        println(date)
    end
    datetime = Date(parse(Int16, date[1:4]), parse(Int8, date[5:6]), parse(Int8, date[7:8]))

    return Dates.value(datetime)
end

"""
    expandTruckMatrices!(nbplannedtrucks, [...])

Return extended matrix of data of planned trucks + extra trucks based on given planned 
trucks data (`_P` matrices).

# Arguments
- `TE`: nbtrucks * nbsuppliers matrix of truck-supplier loading order.
- `TL`: Lengths of each truck.
- `TW`: Width of each truck.
- `TH`: Height of each truck.
- `TKE`: nbtrucks * nbsupplierdocks matrix of truck-supplier dock loading order.
- `TGE`: nbtrucks * nbplants Truck-plant dock loading order matrix.
- `TDA`: Arrival time of each truck.
- `TDE`: ? TODO
- `TU`: nbtrucks * nbsuppliers matrix of candidate suppliers picked-up by each truck.
- `TP`: nbtrucks * nbplants matrix of plant of each truck.
- `TG`: Plant docks of each truck.
- `TK`: nbtrucks * nbsupplierdocks matrix of supplier docks delivered by each truck.
- `TID`: ? TODO
- `TR`: nbtrucks * nbitems matrix of compatible items for each truck.
- `TE_P`: nbplannedtrucks * nbsuppliers matrix of planned truck-supplier loading order.
- `TL_P`: Lengths of each planned truck.
- `TW_P`: Width of each planned truck.
- `TH_P`: Height of each planned truck.
- `TKE_P`: nbplannedtrucks * nbsupplierdocks matrix of planned truck-supplier dock loading order.
- `TGE_P`: nbplannedtrucks * nbplants Truck-plant dock loading order matrix.
- `TDA_P`: Arrival time of each planned truck.
- `TDE_P`: ? TODO
- `TU_P`: nbplannedtrucks * nbsuppliers matrix of candidate suppliers picked-up by each planned truck.
- `TP_P`: nbplannedtrucks * nbplants matrix of plant of each planned truck.
- `TG_P`: Plant docks of each planned truck.
- `TK_P`: nbplannedtrucks * nbsupplierdocks matrix of supplier docks delivered by each planned truck.
- `TR_P`: nbplannedtrucks * nbitems matrix of compatible items for each planned truck.
- `reverse_truckdict`: gives mapping of truck index => planned truck original identifier
- `truckindices`: stores for each planned truck the index of the planned trucks
and indices of corresponding extra trucks.

"""
function expandTruckMatrices!(nbplannedtrucks,
     TE,
     TL,
     TW,
     TH,
     TKE,
     TGE,
     TDA,
     TDE,
     TU,
     TP,
     TG,
     TK,
     TID,
     TR,
     TE_P,
     TL_P,
     TW_P,
     TH_P,
     TKE_P,
     TGE_P,
     TDA_P,
     TDE_P,
     TU_P,
     TP_P,
     TG_P,
     TK_P,
     TR_P,
     reverse_truckdict,
     truckindices
    )
    # j: index of extra trucks
    j = nbplannedtrucks
    # For each planned truck
    for p in 1:nbplannedtrucks
        push!(truckindices[p], p)
        # Get planned truck as first lines of the final matrices
        TE[p, :] .= TE_P[p, :]
        TL[p] = TL_P[p]
        TW[p] = TW_P[p]
        TH[p] = TH_P[p]
        TKE[p, :] .= TKE_P[p, :]
        TGE[p, :] .= TGE_P[p, :]
        TDA[p] = TDA_P[p]
        TDE[p] = TDE_P[p]
        TU[p, :] .= TU_P[p, :]
        TP[p, :] .= TP_P[p, :]
        TG[p, :] .= TG_P[p, :]
        TK[p, :] .= TK_P[p, :]
        TID[p] = reverse_truckdict[p]

        TR[p, :] .= TR_P[p, :]
        
        # for each candidate items, add an extra truck
        for e in 1:sum(TR_P[p,:])-1
            j = j + 1
            push!(truckindices[p], j)

            # Fill relevant truck information
            # Extra trucks have the same data than corresponding planned trucks
            # TODO reduce memory footprint using only `_P` matrices
            tail = "_E" * string(e)
            TID[j] = reverse_truckdict[p] * tail
            TE[j, :] .= TE_P[p, :]
            TL[j] = TL_P[p]
            TW[j] = TW_P[p]
            TH[j] = TH_P[p]
            TKE[j, :] .= TKE_P[p, :]
            TGE[j, :] .= TGE_P[p, :]
            TDA[j] = TDA_P[p]
            TDE[j] = TDE_P[p]
            TU[j, :] .= TU_P[p, :]
            TP[j, :] .= TP_P[p, :]
            TG[j, :] .= TG_P[p, :]
            TK[j, :] .= TK_P[p, :]
            
            TR[j, :] .= TR_P[p, :]
        end
        # @debug begin
        #     sleep(10)
        # end

    end
    # @debug j j 

    

end

"""
    fillItems!([...])

Given matrices of item data, fill these given item data file `"input_items.csv"`
found in `instancepath` folder.

# Arguments

- `IU`: nbitems * nbsuppliers matrix of supplier of each item.
- `IP`: nbitems * nbplants matrix of plant of each item.
- `IK`: nbitems * nbsupplierdocks matrix of supplier docks of each item.
- `IPD`: nbitems * nbplantdocks matrix of plant docks of each item.
- `IDL`: latest arrival time of each item.
- `IDE`: earliest arrival time of each item.
- `IS`: stackability code index of each item.
- `_IO`: forced orientation of each item.
- `IL`: length of each item.
- `IW`: width of each item.
- `IH`: height of each item.
- `stackabilitycodedict`: assigns an index to each stackability code
- `supplierdict`: assigns an index to each supplier's id.
- `plantdict`: assigns an index to each plant's id.
- `supplierdockdict`: assigns an index to each supplier dock's id.
- `plantdockdict`: assigns an index to each plant dock's id.
- `itemdict`: assigns an index to each item's id.
- `instancepath`: path to the folder containing the item data file.
"""
function fillItems!(itemmatrices, 
    stackabilitycodedict, supplierdict, plantdict, supplierdockdict, 
    plantdockdict, itemdict, instancepath, item_copies, inventory_cost, pkg_code,
    max_stackability, stackability_code
)
    nbstackabilitycodes = 0
    open(*(instancepath, "input_items.csv")) do input_itemsfile
        # fill lines of data for each data type and each item.
        for (i, row) in enumerate(CSV.File(input_itemsfile, normalizenames=true, delim=';', decimal=',', stripwhitespace=true, types=String))
            itemmatrices["IU"][i, supplierdict[row[:Supplier_code]]] = 1.0
            itemmatrices["IP"][i, plantdict[row[:Plant_code]]] = 1.0
            # sometimes supplier and plant docks are not given which is not a bug
            itemmatrices["IK"][i, supplierdockdict[row[:Supplier_code]*"__"*(ismissing(row[:Supplier_dock]) ? "missing" : row[:Supplier_dock])]] = 1.0
            itemmatrices["IPD"][i, plantdockdict[row[:Plant_code]*"__"*(ismissing(row[:Plant_dock]) ? "missing" : row[:Plant_dock])]] = 1.0 # TODO since plant docks are unique, no need of plants?
            # TODO some plants don't have docks though
            itemmatrices["IDL"][i] = parse(Float64, row[:Latest_arrival_time])
            itemmatrices["IDE"][i] = parse(Float64, row[:Earliest_arrival_time])
            # No given orientation means no constraint
            itemmatrices["_IO"][i] = row[:Forced_orientation] == "none" ? missing : row[:Forced_orientation] == "widthwise" ? 1 : 0
            # we assign an index to each unique stackability code
            if !haskey(stackabilitycodedict, row[:Stackability_code])
                nbstackabilitycodes = nbstackabilitycodes + 1.0
                stackabilitycodedict[row[:Stackability_code]] = nbstackabilitycodes
            end
            # we use the index of stackability codes in the problem data
            itemmatrices["IS"][i] = stackabilitycodedict[row[:Stackability_code]]
            itemdict[i] = row[:Item_ident]
            itemmatrices["IL"][i] = parse(Float64, replace(row[:Length], "," => "."))
            itemmatrices["IW"][i] = parse(Float64, replace(row[:Width], "," => "."))
            itemmatrices["IH"][i] = parse(Float64, replace(row[:Height], "," => "."))
            itemmatrices["InH"][i] = parse(Float64, replace(row[:Nesting_height], "," => "."))
            itemmatrices["IM"][i] = parse(Float64, replace(row[:Weight], "," => "."))
            item_copies[i] = parse(Float64, row[:Number_of_items])
            pkg_code[i] = row[:Package_code]
            inventory_cost[i] = parse(Float64, row[:Inventory_cost])
            max_stackability[i] = parse(Float64, row[:Max_stackability])
            stackability_code[i] = row[:Stackability_code]

        end
    end
    return nbstackabilitycodes
end

"""
    count!(d::Dict{String, <:Integer}, nb::Integer, rowvalue::String)

Is used to map an index by order of appearance to each `rowvalue` in `d`.

See also [`countTrucks!`](@ref) which makes use of it.
"""
function count!(d::Dict{String, <:Integer}, nb::Integer, rowvalue::String)
    if !haskey(d, rowvalue)
        nb = nb + 1

        d[rowvalue] = nb

    end
    return nb
end

"""
    countTrucks!(truckdict, supplierdict, supplierdockdict, plantdict, plantdockdict, instancepath)

Map an index by order of appearance to each supplier's id, supplier dock id's, 
plant id's and plant dock id's, in given initially empty dictionaries.
Count and return nbplannedtrucks, nbsuppliers, nbsupplierdocks, nbplants, nbplantdocks.

See also [`count!`](@ref).
"""
function countTrucks!(
    truckdict,
    supplierdict,
    supplierdockdict,
    plantdict,
    plantdockdict,
    instancepath)

    nbplannedtrucks = 0
    nbsuppliers = 0
    nbsupplierdocks = 0
    nbplants = 0
    nbplantdocks = 0
    open(*(instancepath, "input_trucks.csv")) do input_trucksfile
        for row in CSV.File(input_trucksfile, normalizenames=true, delim=';', decimal=',', stripwhitespace=true, types=String)
            nbplannedtrucks = count!(truckdict, nbplannedtrucks, row[:Id_truck])
            nbsuppliers = count!(supplierdict, nbsuppliers, row[:Supplier_code])

            # Longer names for docks because dock names can be the same between 2 suppliers or plants
            nbsupplierdocks = count!(supplierdockdict, nbsupplierdocks, row[:Supplier_code]*"__"*(ismissing(row[:Supplier_dock]) ? "missing" : row[:Supplier_dock]))

            nbplants = count!(plantdict, nbplants, row[:Plant_code])

            nbplantdocks = count!(plantdockdict, nbplantdocks, row[:Plant_code]*"__"*(ismissing(row[:Plant_dock]) ? "missing" : row[:Plant_dock]))
            
        end
    end
    @debug nbplannedtrucks nbplannedtrucks
    return nbplannedtrucks, nbsuppliers, nbsupplierdocks, nbplants, nbplantdocks
end

"""
    fillPlannedTruckMatrices!([...])

Given data folder `instancepath`, fill data matrices for each planned truck.

# Arguments

- `TE_P`: nbplannedtrucks * nbsuppliers matrix of planned truck-supplier loading order.
- `TL_P`: Lengths of each planned truck.
- `TW_P`: Width of each planned truck.
- `TH_P`: Height of each planned truck.
- `TKE_P`: nbplannedtrucks * nbsupplierdocks matrix of planned truck-supplier dock loading order.
- `TGE_P`: nbplannedtrucks * nbplants Truck-plant dock loading order matrix.
- `TDA_P`: Arrival time of each planned truck.
- `TDE_P`: ? TODO
- `TU_P`: nbplannedtrucks * nbsuppliers matrix of candidate suppliers picked-up by each planned truck.
- `TP_P`: nbplannedtrucks * nbplants matrix of plant of each planned truck.
- `TG_P`: Plant docks of each planned truck.
- `TK_P`: nbplannedtrucks * nbsupplierdocks matrix of supplier docks delivered by each planned truck.
- `TR_P`: nbplannedtrucks * nbitems matrix of compatible items for each planned truck.
- `instancepath`: path to the folder containing ` "input_trucks.csv"`.
- `nbitems`: number of items.
- `truckdict`: assigns an index to each truck's id.
- `supplierdict`: assigns an index to each supplier's id.
- `supplierdockdict`: assigns an index to each supplier dock's id.
- `plantdict`: assigns an index to each plant's id.
- `plantdockdict`: assigns an index to each plant dock's id.
- `item_productcodes`: stores for each item's index its product code.
"""
function fillPlannedTruckMatrices!(truckmatrices_P, instancepath, nbitems, truckdict, supplierdict, 
    supplierdockdict, plantdict, plantdockdict, item_productcodes, max_stack_weights)

    # For each line of input_trucks:
    open(*(instancepath, "input_trucks.csv")) do input_trucksfile
        for row in CSV.File(input_trucksfile, normalizenames=true, delim=';', decimal=',', stripwhitespace=true, types=String)
            truck_ind = truckdict[row[:Id_truck]]
            # Fill relevant truck information
            truckmatrices_P["TE_P"][
                truck_ind, 
                supplierdict[row[:Supplier_code]]
            ] = parse(Float64, row[:Supplier_loading_order])

            truckmatrices_P["TL_P"][truck_ind] = parse(Float64, row[:Length])
            truckmatrices_P["TW_P"][truck_ind] = parse(Float64, row[:Width])
            truckmatrices_P["TH_P"][truck_ind] = parse(Float64, row[:Height])

            custom_supplierdock_code = 
                row[:Supplier_code]*"__"*(ismissing(row[:Supplier_dock]) ? 
                    "missing" : row[:Supplier_dock])

            truckmatrices_P["TKE_P"][
                truck_ind, supplierdockdict[custom_supplierdock_code]
            ] = parse(Float64, row[:Supplier_dock_loading_order])

            custom_plantdock_code = 
                row[:Plant_code]*"__"*(ismissing(row[:Plant_dock]) ? 
                    "missing" : row[:Plant_dock])

            truckmatrices_P["TGE_P"][
                truck_ind, plantdockdict[custom_plantdock_code]
            ] = parse(Float64, row[:Plant_dock_loading_order])

            truckmatrices_P["TDA_P"][truck_ind] = convert(Float64, to_days(string(parse(Int64, row[:Arrival_time]))))
            truckmatrices_P["TDE_P"][truck_ind] = convert(Float64, to_days(string(parse(Int64, row[:Arrival_time])))) # ? TODO

            truckmatrices_P["TU_P"][truck_ind, supplierdict[row[:Supplier_code]]] = 1.0
            truckmatrices_P["TP_P"][truck_ind, plantdict[row[:Plant_code]]] = 1.0

            truckmatrices_P["TG_P"][
                truck_ind, 
                plantdockdict[row[:Plant_code]*"__"*(ismissing(row[:Plant_dock]) ? 
                    "missing" : row[:Plant_dock])]
            ] = 1.0

            truckmatrices_P["TK_P"][
                truck_ind, supplierdockdict[row[:Supplier_code]*"__"*(ismissing(row[:Supplier_dock]) ?
                    "missing" : row[:Supplier_dock])]
            ] = 1.0
            
            truckmatrices_P["TEM_P"][truck_ind] = parse(Float64, replace(row[:Max_density], "," => "."))
            
            if !haskey(max_stack_weights, row[:Id_truck])
                max_stack_weights[row[:Id_truck]] = Dict()
            end

            product_code = row[:Product_code]
            max_stack_weights[row[:Id_truck]][product_code] = parse(Float64, replace(row[:Max_weight_on_the_bottom_item_in_stacks], "," => "."))
            
            truckmatrices_P["TMm_P"][truck_ind] = parse(Float64, replace(row[:Max_weight], "," => "."))

            truckmatrices_P["CM_P"][truck_ind] = parse(Float64, replace(row[:CM], "," => "."))
            truckmatrices_P["CJfm_P"][truck_ind] = parse(Float64, replace(row[:CJfm], "," => "."))
            truckmatrices_P["CJfc_P"][truck_ind] = parse(Float64, replace(row[:CJfc], "," => "."))
            truckmatrices_P["CJfh_P"][truck_ind] = parse(Float64, replace(row[:CJfh], "," => "."))
            truckmatrices_P["EM_P"][truck_ind] = parse(Float64, replace(row[:EM], "," => "."))
            truckmatrices_P["EJhr_P"][truck_ind] = parse(Float64, replace(row[:EJhr], "," => "."))
            truckmatrices_P["EJcr_P"][truck_ind] = parse(Float64, replace(row[:EJcr], "," => "."))
            truckmatrices_P["EJeh_P"][truck_ind] = parse(Float64, replace(row[:EJeh], "," => "."))
            truckmatrices_P["EMmr_P"][truck_ind] = parse(Float64, replace(row[:EMmr], "," => "."))
            truckmatrices_P["EMmm_P"][truck_ind] = parse(Float64, replace(row[:EMmm], "," => "."))
            truckmatrices_P["Cost_P"][truck_ind] = parse(Float64, replace(row[:Cost], "," => "."))

            # For each line of input trucks, retrieve Product code. For all items, get 
            # indices of items of same product code, and use it to fill TR
            for i in 1:nbitems
                #Takes a lot of times because trucks * items
                truckmatrices_P["TR_P"][truck_ind, i] = item_productcodes[i] == product_code ? 1.0 : truckmatrices_P["TR_P"][truckdict[row[:Id_truck]], i]
            end
            # @debug begin
            #     if sum([item_productcodes[i] == product_code ? 1.0 : 0.0 for i in 1:nbitems]) >= 1
            #         @debug product_code product_code
            #         @debug "[item_productcodes[i] == product_code ? 1.0 : 0.0 for i in 1:nbitems]" [item_productcodes[i] == product_code ? 1.0 : 0.0 for i in 1:10]
            #         @debug "TR_P[truckdict[row[:Id_truck]], :]" TR_P[truckdict[row[:Id_truck]], 1:10]
            #         sleep(20)
            #     end
            # end
        end
    end
end

"""
    loadinstance(instancepath)

Return all relevant data of truck stack item assignment problem defined in data files found in folder 
`instancepath`.
"""
function loadinstance(instancepath; onlyplanned=false)

    """
    Important notes
    
    The order of items is the one of input_items.
    """
    ## Count items and associate product code to each item's index
    nbitems = 0
    item_productcodes = Vector{String}()
    open(*(instancepath, "input_items.csv")) do input_itemsfile
        for row in CSV.File(input_itemsfile, normalizenames=true, delim=';', decimal=',', stripwhitespace=true, types=String)
            nbitems = nbitems + 1
            # @debug "row[:Product_code]" row[:Product_code]
            # @debug "row[:Product_code] type" typeof(row[:Product_code])
            push!(item_productcodes, String(row[:Product_code]))
        end
    end
    
    """
    The truck order will be the one of the appearance of each truck in input_trucks
    """
    
    # Make list of trucks with corresponding dictionary
    truckdict = Dict{String, Int64}()

    # Make list of suppliers with corresponding dictionary
    supplierdict = Dict{String, Int64}()

    # Make list of supplier docks with corresponding dictionary
    supplierdockdict = Dict{String, Int64}()

    # Make list of plants with corresponding dictionary
    plantdict = Dict{String, Int64}()

    # Make list of plant docks with corresponding dictionary
    plantdockdict = Dict{String, Int64}()
    
    
    ## Count planned trucks, suppliers, supplier docks, plants, plant docks
    ## Make dictionary allowing to find a truck/supplier/supplier dock/plant/plant 
    ## dock's index from its  string ID
    nbplannedtrucks, nbsuppliers, nbsupplierdocks, nbplants, nbplantdocks  = countTrucks!(truckdict, supplierdict, 
    supplierdockdict, plantdict, plantdockdict, instancepath)    
    
    TE_P = Matrix{Union{Float64, Missing}}(missing, nbplannedtrucks, nbsuppliers)
    
    TL_P = Vector{Union{Float64, Missing}}(missing, nbplannedtrucks)
    TW_P = Vector{Union{Float64, Missing}}(missing, nbplannedtrucks)
    TH_P = Vector{Union{Float64, Missing}}(missing, nbplannedtrucks)
    
    TKE_P = Matrix{Union{Float64, Missing}}(missing, nbplannedtrucks, nbsupplierdocks)
    fill!(TKE_P, nbsupplierdocks)
    
    TGE_P = Matrix{Union{Float64, Missing}}(missing, nbplannedtrucks, nbplantdocks)
    fill!(TGE_P, nbplantdocks)
    
    TDA_P = Vector{Union{Float64, Missing}}(missing, nbplannedtrucks)
    TDE_P = Vector{Union{Float64, Missing}}(missing, nbplannedtrucks)
    
    TU_P = falses(nbplannedtrucks, nbsuppliers)
    TP_P = falses(nbplannedtrucks, nbplants)
    TK_P = falses(nbplannedtrucks, nbsupplierdocks)
    
    TG_P = falses(nbplannedtrucks, nbplantdocks)
    
    
    TR_P = falses(nbplannedtrucks, nbitems) # TR is expanded, it will contain also only items which docks are delivered by the truck
    
    TMm_P = Vector{Union{Float64, Missing}}(missing, nbplannedtrucks)
    TEM_P = Vector{Union{Float64, Missing}}(missing, nbplannedtrucks)
    max_stack_weights_P = Dict{String, Any}()

    EMmm_P = Vector{Union{Float64, Missing}}(missing, nbplannedtrucks)
    EMmr_P = Vector{Union{Float64, Missing}}(missing, nbplannedtrucks)
    CM_P = Vector{Union{Float64, Missing}}(missing, nbplannedtrucks)
    CJfm_P = Vector{Union{Float64, Missing}}(missing, nbplannedtrucks)
    CJfc_P = Vector{Union{Float64, Missing}}(missing, nbplannedtrucks)
    CJfh_P = Vector{Union{Float64, Missing}}(missing, nbplannedtrucks)
    EM_P = Vector{Union{Float64, Missing}}(missing, nbplannedtrucks)
    EJhr_P = Vector{Union{Float64, Missing}}(missing, nbplannedtrucks)
    EJcr_P = Vector{Union{Float64, Missing}}(missing, nbplannedtrucks)
    EJeh_P = Vector{Union{Float64, Missing}}(missing, nbplannedtrucks)
    Cost_P = Vector{Union{Float64, Missing}}(missing, nbplannedtrucks)
    truckmatrices_P = Dict()
    
    truckmatrices_P["TE_P"] = TE_P
    truckmatrices_P["TL_P"] = TL_P
    truckmatrices_P["TW_P"] = TW_P
    truckmatrices_P["TH_P"] = TH_P
    truckmatrices_P["TKE_P"] = TKE_P
    truckmatrices_P["TGE_P"] = TGE_P
    truckmatrices_P["TDA_P"] = TDA_P
    truckmatrices_P["TDE_P"] = TDE_P
    truckmatrices_P["TU_P"] = TU_P
    truckmatrices_P["TP_P"] = TP_P
    truckmatrices_P["TK_P"] = TK_P
    truckmatrices_P["TG_P"] = TG_P
    truckmatrices_P["TR_P"] = TR_P
    truckmatrices_P["TMm_P"] = TMm_P
    truckmatrices_P["TEM_P"] = TEM_P
    truckmatrices_P["max_stack_weights_P"] = max_stack_weights_P
    truckmatrices_P["EMmm_P"] = EMmm_P
    truckmatrices_P["EMmr_P"] = EMmr_P
    truckmatrices_P["CM_P"] = CM_P
    truckmatrices_P["CJfm_P"] = CJfm_P
    truckmatrices_P["CJfc_P"] = CJfc_P
    truckmatrices_P["CJfh_P"] = CJfh_P
    truckmatrices_P["EM_P"] = EM_P
    truckmatrices_P["EJhr_P"] = EJhr_P
    truckmatrices_P["EJcr_P"] = EJcr_P
    truckmatrices_P["EJeh_P"] = EJeh_P
    truckmatrices_P["Cost_P"] = Cost_P

    ## Fill planned truck matrices with info in input_trucksfile
    fillPlannedTruckMatrices!(truckmatrices_P, instancepath, nbitems, truckdict, 
    supplierdict, supplierdockdict, plantdict, plantdockdict, item_productcodes, max_stack_weights_P)
    
    IU = falses(nbitems, nbsuppliers)
    IP = falses(nbitems, nbplants)
    IK = falses(nbitems, nbsupplierdocks)
    IPD = falses(nbitems, nbplantdocks)
    
    _IO = Vector{Union{Float64, Missing}}(missing, nbitems)

    IL = Vector{Float64}(undef, nbitems)
    IW = Vector{Float64}(undef, nbitems)
    IH = Vector{Float64}(undef, nbitems)
    InH = Vector{Float64}(undef, nbitems)
    IM = Vector{Float64}(undef, nbitems)
    IS = Vector{Union{Float64, Missing}}(missing, nbitems)
    pkg_code = Vector{Union{String, Missing}}(missing, nbitems)
    IDL = Vector{Union{Float64, Missing}}(missing, nbitems)
    stackability_code = Vector{Union{String, Missing}}(missing, nbitems)
    IDE = Vector{Union{Float64, Missing}}(missing, nbitems)
    itemdict = Dict{Int64, String}()
    item_copies = Vector{Union{Float64, Missing}}(missing, nbitems)
    inventory_cost = Vector{Union{Float64, Missing}}(missing, nbitems)
    max_stackability = Vector{Union{Float64, Missing}}(missing, nbitems)
    ## Fill Item info matrices with data from input_items file
    stackabilitycodedict = Dict{String, Float64}()

    itemmatrices = Dict()
    itemmatrices["IU"] = IU
    itemmatrices["IP"] = IP
    itemmatrices["IK"] = IK
    itemmatrices["IPD"] = IPD
    itemmatrices["_IO"] = _IO
    itemmatrices["IL"] = IL
    itemmatrices["IW"] = IW
    itemmatrices["IH"] = IH
    itemmatrices["InH"] = InH
    itemmatrices["IM"] = IM
    itemmatrices["IS"] = IS
    itemmatrices["IDL"] = IDL
    itemmatrices["IDE"] = IDE

    nbstackabilitycodes = fillItems!(itemmatrices, stackabilitycodedict, 
        supplierdict, plantdict, supplierdockdict, 
        plantdockdict, itemdict, instancepath, item_copies, inventory_cost, pkg_code,
        max_stackability, stackability_code
    )
    
    # Expand TR with information about docks
    # For each truck, for each item, if the truck doesn't stop at the supplier & supplier dock of the item or 
    # if it doesn't stop by the plant & plant dock of the item: replace with 0
    @debug "sum(TR_P)" sum(TR_P)
    
    ## Adjust TR with additional information
    # This loop might be useless or too restrictive, because it seems that the candidate list of 
    # each truck already includes this type of info TODO
    for t in 1:nbplannedtrucks
        for i in 1:nbitems
            # if item `i` is compatible with truck `t` according to the candidate list
            if TR_P[t,i] == 1
                # if the supplier dock of i is not in the supplier docks of t
                # or the same but with plant docks
                # or the truck arrives after the latest arrival date allowed for i 
                # or the truck's plant is not the item's
                # or the item's supplier is not in the truck's
                if !*((IK[i,:] .<= TK_P[t,:])...) || 
                    !*((IPD[i,:] .<= TG_P[t,:])...) || 
                    TDA_P[t] > IDL[i] ||
                    TP_P[t] != IP[i] ||
                    !*((IU[i,:] .<= TU_P[t,:])...)
                    # Remove i from compatible items for t
                    TR_P[t,i] = 0.0
                end
            end
        end
    end
    
    # @debug "" !*((IK[5,:] .<= TK_P[truckdict["P192711301"],:])...) || !*((IPD[5,:] .<= TG_P[truckdict["P192711301"],:])...) || TDA_P[truckdict["P192711301"]] > IDL[5]
    # @debug begin
        #     println(!*((IK[5,:] .<= TK_P[truckdict["P192711301"],:])...))
        #     println(!*((IPD[5,:] .<= TG_P[truckdict["P192711301"],:])...))
        #     println(TDA_P[truckdict["P192711301"]] > IDL[5])
        #     println("Item ident of 5", itemdict[5])
    # end
    TE = nothing
    TL = nothing
    TW = nothing
    TH = nothing
    TKE = nothing
    TGE = nothing
    TDA = nothing
    TDE = nothing
    TU = nothing
    TP = nothing
    TK = nothing
    TG = nothing
    TR = nothing
    TID = nothing
    nbtrucks = nbplannedtrucks
    reverse_truckdict = nothing
    truckindices = nothing
    reverse_truckdict = Dict(value => key for (key, value) in truckdict)
    if !onlyplanned
        ## Compute actual number of trucks (planned trucks + extra ones)
        # The total number of trucks could be nbplannedtrucks * nbitems, but a 
        # smarter way would be to have
        # nbtrucks = sum(TR) # sum of number of candidate items for each truck
        # nbtrucks = nbplannedtrucks * nbitems # Very bad idea (millions of trucks)
        nbtrucks = nbplannedtrucks + sum([sum(TR_P[t, :]) > 0 ? sum(TR_P[t, :])-1 : 0 for t in 1:nbplannedtrucks] )
        # @debug nbplannedtrucks nbplannedtrucks
        # @debug sum(TR_P) sum(TR_P)
        # @debug nbitems nbitems
        # @debug "" nbtrucks
        # @debug "nbplannedtrucks * nbitems" nbplannedtrucks * nbitems
        
        
        TE = Matrix{Union{Float64, Missing}}(missing, nbtrucks, nbsuppliers)
        
        TL = Vector{Union{Float64, Missing}}(missing, nbtrucks)
        TW = Vector{Union{Float64, Missing}}(missing, nbtrucks)
        TH = Vector{Union{Float64, Missing}}(missing, nbtrucks)
        
        TKE = Matrix{Union{Float64, Missing}}(missing, nbtrucks, nbsupplierdocks)
        fill!(TKE, nbsupplierdocks)
        
        TGE = Matrix{Union{Float64, Missing}}(missing, nbtrucks, nbplantdocks)
        fill!(TGE, nbplantdocks)
        
        TDA = Vector{Union{Float64, Missing}}(missing, nbtrucks)
        TDE = Vector{Union{Float64, Missing}}(missing, nbtrucks)
        
        TU = falses(nbtrucks, nbsuppliers)
        TP = falses(nbtrucks, nbplants)
        TK = falses(nbtrucks, nbsupplierdocks)
        
        TG = falses(nbtrucks, nbplantdocks)
        
        
        TR = falses(nbtrucks, nbitems) # TR is expanded, it will contain also only items which docks are delivered by the truck
        TID = Vector{Union{String, Missing}}(missing, nbtrucks)
        
        # Associate to each planned truck, the indices of the planned truck + extra trucks
        truckindices = [Vector{Integer}() for p in 1:nbplannedtrucks]

        ## Fill the actual truck matrices from the planned trucks (expand the planned 
        ## trucks matrix with extra trucks)
        expandTruckMatrices!(nbplannedtrucks, TE, TL, TW, TH, TKE, TGE, TDA, TDE, TU, TP, TG, TK, TID, TR, 
        TE_P, TL_P, TW_P, TH_P, TKE_P, TGE_P, TDA_P, TDE_P, TU_P, TP_P, TG_P, TK_P, TR_P, 
        reverse_truckdict, truckindices)
    end

    coefcostinventory = 0.0
    costtransportation = 0.0
    coefcostextratruck = 0.0
    timelimit = 0.0
    open(*(instancepath, "input_parameters.csv")) do file
        for row in CSV.File(file, normalizenames=true, delim=';', decimal=',', stripwhitespace=true, types=String)
            coefcostinventory = parse(Float64, replace(row[:Coefficient_inventory_cost], "," => "."))
            costtransportation = parse(Float64, replace(row[:Coefficient_transportation_cost], "," => "."))
            coefcostextratruck = parse(Float64, replace(row[:Coefficient_cost_extra_truck], "," => "."))
            timelimit = row[CSV.normalizename("timelimit_(sec)")]       

        end
    end
    result = Dict{Any, Any}()
    result["item_productcodes"] = item_productcodes
    result["truckdict"] = truckdict
    result["supplierdict"] = supplierdict
    result["supplierdockdict"] = supplierdockdict
    result["plantdict"] = plantdict
    result["plantdockdict"] = plantdockdict
    result["nbplannedtrucks"] = nbplannedtrucks
    result["nbitems"] = nbitems
    result["nbsuppliers"] = nbsuppliers
    result["nbsupplierdocks"] = nbsupplierdocks
    result["nbplants"] = nbplants
    result["nbplantdocks"] = nbplantdocks
    result["TE_P"] = TE_P
    result["TL_P"] = TL_P
    result["TW_P"] = TW_P
    result["TH_P"] = TH_P
    result["TKE_P"] = TKE_P
    result["TGE_P"] = TGE_P
    result["TDA_P"] = TDA_P
    result["TDE_P"] = TDE_P
    result["TU_P"] = TU_P
    result["TP_P"] = TP_P
    result["TK_P"] = TK_P
    result["TG_P"] = TG_P
    result["TR_P"] = TR_P
    result["CM_P"] = CM_P
    result["TEM_P"] = TEM_P
    result["TMm_P"] = TMm_P
    result["CM_P"] = CM_P
    result["CJfm_P"] = CJfm_P
    result["CJfc_P"] = CJfc_P
    result["CJfh_P"] = CJfh_P
    result["EM_P"] = EM_P
    result["EJhr_P"] = EJhr_P
    result["EJcr_P"] = EJcr_P
    result["EJeh_P"] = EJeh_P
    result["EMmr_P"] = EMmr_P
    result["EMmm_P"] = EMmm_P
    result["Cost_P"] = Cost_P
    result["max_stack_weights_P"] = max_stack_weights_P

    result["IU"] = IU
    result["IP"] = IP
    result["IK"] = IK
    result["IPD"] = IPD
    result["IS"] = IS
    result["_IO"] = _IO
    result["IL"] = IL
    result["IW"] = IW
    result["IH"] = IH
    result["InH"] = InH
    result["IM"] = IM
    result["IDL"] = IDL
    result["IDE"] = IDE
    result["stackabilitycodedict"] = stackabilitycodedict
    result["nbtrucks"] = nbtrucks
    result["TE"] = TE
    result["TL"] = TL
    result["TW"] = TW
    result["TH"] = TH
    result["TKE"] = TKE
    result["TGE"] = TGE
    result["TDA"] = TDA
    result["TDE"] = TDE
    result["TU"] = TU
    result["TP"] = TP
    result["TK"] = TK
    result["TG"] = TG
    result["TR"] = TR
    result["TID"] = TID
    result["reverse_truckdict"] = reverse_truckdict
    result["truckindices"] = truckindices
    result["coefcostinventory"] = coefcostinventory
    result["costtransportation"] = costtransportation
    result["coefcostextratruck"] = coefcostextratruck
    result["timelimit"] = timelimit
    result["item_copies"] = item_copies
    result["nbitems"] = nbitems
    result["pkg_code"] = pkg_code
    result["itemdict"] = itemdict
    result["max_stackability"] = max_stackability
    result["inventory_cost"] = inventory_cost
    result["stackability_code"] = stackability_code

    return result
    
end
# end
