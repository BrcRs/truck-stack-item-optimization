# module TSIInstanceLoader

# using LinearAlgebra
# using JuMP
# using Cbc
using CSV

# using FilePaths
using Logging
# using Profile

# using ArgParse

logger = ConsoleLogger(stdout, Logging.Debug, show_limited=false)

global_logger(logger)





function expandTruckMatrices!(nbplannedtrucks,
     TE,
     TL,
     TW,
     TH,
     TKE,
     TGE,
     TDA,
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
     TU_P,
     TP_P,
     TG_P,
     TK_P,
     TR_P,
     reverse_truckdict
    )

    # For each planned truck
    j = nbplannedtrucks
    for p in 1:nbplannedtrucks

        # Get planned truck as first lines of the bunch
        TE[p, :] .= TE_P[p, :]
        TL[p] = TL_P[p]
        TW[p] = TW_P[p]
        TH[p] = TH_P[p]
        TKE[p, :] .= TKE_P[p, :]
        TGE[p, :] .= TGE_P[p, :]
        TDA[p] = TDA_P[p]
        TU[p, :] .= TU_P[p, :]
        TP[p, :] .= TP_P[p, :]
        TG[p, :] .= TG_P[p, :]
        TK[p, :] .= TK_P[p, :]
        TID[p] = reverse_truckdict[p]

        TR[p, :] .= TR_P[p, :]
        
        # for each candidate items, add an extra truck
        for e in 1:sum(TR_P[p,:])-1
            j = j + 1

            # Fill relevant truck information
            tail = "_E" * string(e)
            TID[j] = reverse_truckdict[p] * tail
            TE[j, :] .= TE_P[p, :]
            TL[j] = TL_P[p]
            TW[j] = TW_P[p]
            TH[j] = TH_P[p]
            TKE[j, :] .= TKE_P[p, :]
            TGE[j, :] .= TGE_P[p, :]
            TDA[j] = TDA_P[p]
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

function fillItems!(IU, IP, IK, IPD, IDL, IDE, IS, _IO, IL, IW, IH, stackabilitycodedict, supplierdict, plantdict, supplierdockdict, plantdockdict, itemdict, instancepath)
    nbstackabilitycodes = 0
    open(*(instancepath, "input_items.csv")) do input_itemsfile

        for (i, row) in enumerate(CSV.File(input_itemsfile, normalizenames=true, delim=';', decimal=',', stripwhitespace=true, types=String))
            IU[i, supplierdict[row[:Supplier_code]]] = 1.0
            IP[i, plantdict[row[:Plant_code]]] = 1.0
            IK[i, supplierdockdict[row[:Supplier_code]*"__"*(ismissing(row[:Supplier_dock]) ? "missing" : row[:Supplier_dock])]] = 1.0
            IPD[i, plantdockdict[row[:Plant_code]*"__"*(ismissing(row[:Plant_dock]) ? "missing" : row[:Plant_dock])]] = 1.0
            IDL[i] = parse(Float64, row[:Latest_arrival_time])
            IDE[i] = parse(Float64, row[:Earliest_arrival_time])
            _IO[i] = row[:Forced_orientation] == "none" ? missing : row[:Forced_orientation] == "widthwise" ? 1 : 0
            if !haskey(stackabilitycodedict, row[:Stackability_code])
                nbstackabilitycodes = nbstackabilitycodes + 1.0
                stackabilitycodedict[row[:Stackability_code]] = nbstackabilitycodes
            end
            IS[i] = stackabilitycodedict[row[:Stackability_code]]
            itemdict[i] = row[:Item_ident]
            IL[i] = parse(Float64, row[:Length])
            IW[i] = parse(Float64, row[:Width])
            IH[i] = parse(Float64, row[:Height])
        end
    end
    return nbstackabilitycodes
end

function count!(d, nb, rowvalue)
    if !haskey(d, rowvalue)
        nb = nb + 1

        d[rowvalue] = nb

    end
    return nb
end

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

function fillPlannedTruckMatrices!(TE_P, TL_P, TW_P, TH_P, TKE_P, TGE_P, TDA_P, 
    TU_P, TP_P, TG_P, TK_P, TR_P, instancepath, nbitems, truckdict, supplierdict, 
    supplierdockdict, plantdict, plantdockdict, item_productcodes)
    # For each line of input_trucks:
    open(*(instancepath, "input_trucks.csv")) do input_trucksfile
        for row in CSV.File(input_trucksfile, normalizenames=true, delim=';', decimal=',', stripwhitespace=true, types=String)
            
            # Fill relevant truck information
            TE_P[truckdict[row[:Id_truck]], supplierdict[row[:Supplier_code]]] = parse(Float64, row[:Supplier_loading_order])
            TL_P[truckdict[row[:Id_truck]]] = parse(Float64, row[:Length])
            TW_P[truckdict[row[:Id_truck]]] = parse(Float64, row[:Width])
            TH_P[truckdict[row[:Id_truck]]] = parse(Float64, row[:Height])
            TKE_P[truckdict[row[:Id_truck]], supplierdockdict[row[:Supplier_code]*"__"*(ismissing(row[:Supplier_dock]) ? "missing" : row[:Supplier_dock])]] = parse(Float64, row[:Supplier_dock_loading_order])
            TGE_P[truckdict[row[:Id_truck]], plantdockdict[row[:Plant_code]*"__"*(ismissing(row[:Plant_dock]) ? "missing" : row[:Plant_dock])]] = parse(Float64, row[:Plant_dock_loading_order])
            TDA_P[truckdict[row[:Id_truck]]] = parse(Float64, row[:Arrival_time])
            TU_P[truckdict[row[:Id_truck]], supplierdict[row[:Supplier_code]]] = 1.0
            TP_P[truckdict[row[:Id_truck]], plantdict[row[:Plant_code]]] = 1.0


            TG_P[truckdict[row[:Id_truck]], plantdockdict[row[:Plant_code]*"__"*(ismissing(row[:Plant_dock]) ? "missing" : row[:Plant_dock])]] = 1.0
            TK_P[truckdict[row[:Id_truck]], supplierdockdict[row[:Supplier_code]*"__"*(ismissing(row[:Supplier_dock]) ? "missing" : row[:Supplier_dock])]] = 1.0
            
            
            # For each line of input trucks, retrieve Product code. For all items, get 
            # indices of item of same product code, and use it to fill TR
            product_code = row[:Product_code]
            for i in 1:nbitems
                TR_P[truckdict[row[:Id_truck]], i] = item_productcodes[i] == product_code ? 1.0 : TR_P[truckdict[row[:Id_truck]], i]
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

function loadinstance(instancepath)

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
    
    TU_P = falses(nbplannedtrucks, nbsuppliers)
    TP_P = falses(nbplannedtrucks, nbplants)
    TK_P = falses(nbplannedtrucks, nbsupplierdocks)
    
    TG_P = falses(nbplannedtrucks, nbplantdocks)
    
    
    TR_P = falses(nbplannedtrucks, nbitems) # TR is expanded, it will contain also only items which docks are stopped by by the truck
    
    ## Fill planned truck matrices with info in input_trucksfile
    fillPlannedTruckMatrices!(TE_P, TL_P, TW_P, TH_P, TKE_P, TGE_P, TDA_P, 
    TU_P, TP_P, TG_P, TK_P, TR_P, instancepath, nbitems, truckdict, supplierdict, 
    supplierdockdict, plantdict, plantdockdict, item_productcodes)
    
    IU = falses(nbitems, nbsuppliers)
    IP = falses(nbitems, nbplants)
    IK = falses(nbitems, nbsupplierdocks)
    IPD = falses(nbitems, nbplantdocks)
    
    _IO = Vector{Union{Float64, Missing}}(missing, nbitems)

    IL = Vector{Float64}(undef, nbitems)
    IW = Vector{Float64}(undef, nbitems)
    IH = Vector{Float64}(undef, nbitems)
    IS = Vector{Union{Float64, Missing}}(missing, nbitems)
    
    IDL = Vector{Union{Float64, Missing}}(missing, nbitems)
    
    IDE = Vector{Union{Float64, Missing}}(missing, nbitems)
    itemdict = Dict{Int64, String}()
    
    ## Fill Item info matrices with data from input_items file
    stackabilitycodedict = Dict{String, Float64}()
    nbstackabilitycodes = fillItems!(IU, IP, IK, IPD, IDL, IDE, IS, _IO, IL, IW, IH, stackabilitycodedict, supplierdict, plantdict, supplierdockdict, plantdockdict, itemdict, instancepath)
    
    # Expand TR with information about docks
    # For each truck, for each item, if the truck doesn't stop at the supplier & supplier dock of the item or 
    # it doesn't stop by the plant & plant dock of the item: replace with 0
    @debug sum(TR_P) sum(TR_P)
    
    ## Adjust TR with additional information
    # This loop might be useless or too restrictive, because it seems that the candidate list of 
    # each truck already includes this type of info 
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
    
    TU = falses(nbtrucks, nbsuppliers)
    TP = falses(nbtrucks, nbplants)
    TK = falses(nbtrucks, nbsupplierdocks)
    
    TG = falses(nbtrucks, nbplantdocks)
    
    
    TR = falses(nbtrucks, nbitems) # TR is expanded, it will contain also only items which docks are stopped by by the truck
    TID = Vector{Union{String, Missing}}(missing, nbtrucks)
    reverse_truckdict = Dict(value => key for (key, value) in truckdict)

    ## Fill the actual truck matrices from the planned trucks (expand the planned 
    ## trucks matrix with extra trucks)
    expandTruckMatrices!(nbplannedtrucks, TE, TL, TW, TH, TKE, TGE, TDA, TU, TP, TG, TK, TID, TR, 
    TE_P, TL_P, TW_P, TH_P, TKE_P, TGE_P, TDA_P, TU_P, TP_P, TG_P, TK_P, TR_P, 
    reverse_truckdict)
    
    costinventory = 0.0
    costtransportation = 0.0
    costextratruck = 0.0
    timelimit = 0.0
    open(*(instancepath, "input_parameters.csv")) do file
        for row in CSV.File(file, normalizenames=true, delim=';', decimal=',', stripwhitespace=true, types=String)
            costinventory = parse(Float64, replace(row[:Coefficient_inventory_cost], "," => "."))
            costtransportation = parse(Float64, replace(row[:Coefficient_transportation_cost], "," => "."))
            costextratruck = parse(Float64, replace(row[:Coefficient_cost_extra_truck], "," => "."))
            timelimit = row[CSV.normalizename("timelimit_(sec)")]       

        end
    end
    
    return item_productcodes, truckdict, supplierdict, supplierdockdict, plantdict, 
    plantdockdict, nbplannedtrucks, nbitems, nbsuppliers, nbsupplierdocks, nbplants, nbplantdocks,
    TE_P, TL_P, TW_P, TH_P, TKE_P, TGE_P, TDA_P, TU_P, TP_P, TK_P, TG_P, TR_P,
    IU, IP, IK, IPD, IS, _IO, IL, IW, IH, IDL, IDE, stackabilitycodedict, nbtrucks, TE, TL, TW, TH,
    TKE, TGE, TDA, TU, TP, TK, TG, TR, TID, reverse_truckdict, 
    costinventory, costtransportation, costextratruck, timelimit
    
end
# end
