module TSIInstanceLoader

using LinearAlgebra
using JuMP
using Cbc
using CSV

using FilePaths
using Logging
using Profile

using ArgParse

logger = ConsoleLogger(stdout, Logging.Debug, show_limited=false)

global_logger(logger)


function fillXi1!(Xi::BitArray)
    n = size(Xi)[2]
    m = 1
    i = 1
    while n > 0
        Xi[m:m+n-1, i] .= 1
        m = m + n
        n = n - 1
        i = i + 1
    end
end

function identityMat(n)
    mat = falses(n, n)
    for i in 1:n
        mat[i,i] = 1
    end

    return mat
end

function fillXi2!(Xi::BitArray)
    n = size(Xi)[2]
    m = 1
    i = 1
    while n > 0
        print(m, ":")
        println(n)
        # Xi[m:m+n-1, i] .= 1
        # Xi[m:m+n-1, i:size(Xi)[2]] = diagm([1 for _ in 1:n])
        # Xi[m:m+n-1, i:size(Xi)[2]] = identityMat(n)
        Xi[m:m+n-1, i:size(Xi)[2]] = I(n)
        m = m + n
        n = n - 1
        i = i + 1
        println()
    end
end

function ones(n::Int)
    return ones(Int8, n, 1)
end


"""
normalizeValues([1, 8, 15, 8]) -> [1, 2, 3, 2]
"""
function normalizeValues(dataraw)

    datadict = Dict{T, Float64}()
    data = Array{Union{Float64, Missing}}(missing, size(dataraw)[1], 1)

    n = 1
    for (i, v) in enumerate(dataraw)
        
        println(i, " ", v)
        if haskey(datadict, v)
            data[i,1] = datadict[v]
        else
            data[i,1] = n
            datadict[v] = n
            n = n + 1
        end
    end
    return data
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        ## Arg Parse examples
        # "--opt1"
        #     help = "an option with an argument"
        # "--opt2", "-o"
        #     help = "another option with an argument"
        #     arg_type = Int
        #     default = 0
        # "--flag1"
        #     help = "an option without argument, i.e. a flag"
        #     action = :store_true
        "instancePath"
            help = "Path to the folder of the instance to solve."
            required = true
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end
    
    # instancePath = "Instances/AS/"
    instancePath = parsed_args["instancePath"]


    # input_itemsCSV
    # input_parametersCSV
    # input_trucksCSV
    
    """
    Important notes
    
    The order of items is the one if input_items.
    """
    # input_itemsfile = open(*(instancePath, "input_items.csv"))    
    # input_parametersfile = open(*(instancePath, "input_parameters.csv"))    
    # input_trucksfile = open(*(instancePath, "input_trucks.csv"))    

    # input_itemsCSV = CSV.File(input_itemsfile, normalizenames=true, delim=';', decimal=',', stripwhitespace=true, types=String)
    # input_parametersCSV = CSV.File(input_parametersfile, normalizenames=true, delim=';', decimal=',', stripwhitespace=true, types=String)
    # input_trucksCSV = CSV.File(input_trucksfile, normalizenames=true, delim=';', decimal=',', stripwhitespace=true, types=String)
    # input_parametersCSV = CSV.read(open(Path(instancePath, "input_parameters.csv")), normalizenames=true, delim=";", decimal=",", stripwhitespace=true, types=String)
    # input_trucksCSV = CSV.read(open(Path(instancePath, "input_trucks.csv")), normalizenames=true, delim=";", decimal=",", stripwhitespace=true, types=String)


    nbItems = 0
    Iproduct_code = Vector{String}()
    open(*(instancePath, "input_items.csv")) do input_itemsfile
        for row in CSV.File(input_itemsfile, normalizenames=true, delim=';', decimal=',', stripwhitespace=true, types=String)
            nbItems = nbItems + 1
            # @debug "row[:Product_code]" row[:Product_code]
            # @debug "row[:Product_code] type" typeof(row[:Product_code])
            push!(Iproduct_code, String(row[:Product_code]))
        end
    end

    """
    The truck order will be the one of the appearance of each truck in input_trucks
    """
    
    # Make list of trucks with corresponding dictionary
    truckDict = Dict{String, Int64}()
    # Make list of suppliers with corresponding dictionary
    supplierDict = Dict{String, Int64}()
    # Make list of supplier docks with corresponding dictionary
    supplierDockDict = Dict{String, Int64}()
    # Make list of plants with corresponding dictionary
    plantDict = Dict{String, Int64}()
    # Make list of plant docks with corresponding dictionary
    plantDockDict = Dict{String, Int64}()
    
    nbPlannedTrucks = 0
    nbSuppliers = 0
    nbSupplierDocks = 0
    nbPlants = 0
    nbPlantDocks = 0
    open(*(instancePath, "input_trucks.csv")) do input_trucksfile
        for row in CSV.File(input_trucksfile, normalizenames=true, delim=';', decimal=',', stripwhitespace=true, types=String)
            if !haskey(truckDict, row[:Id_truck])
                nbPlannedTrucks = nbPlannedTrucks + 1
                truckDict[row[:Id_truck]] = nbPlannedTrucks
            end
            if !haskey(supplierDict, row[:Supplier_code])
                nbSuppliers = nbSuppliers + 1
                # @debug "" typeof(row[:Supplier_code])
                supplierDict[row[:Supplier_code]] = nbSuppliers
            end
            # Longer names for docks because dock names can be the same between 2 suppliers or plants
            if !haskey(supplierDockDict, row[:Supplier_code]*"__"*(ismissing(row[:Supplier_dock]) ? "missing" : row[:Supplier_dock]))
                nbSupplierDocks = nbSupplierDocks + 1

                supplierDockDict[row[:Supplier_code]*"__"*(ismissing(row[:Supplier_dock]) ? "missing" : row[:Supplier_dock])] = nbSupplierDocks

            end
            if !haskey(plantDict, row[:Plant_code])
                nbPlants = nbPlants + 1
                plantDict[row[:Plant_code]] = nbPlants
            end
            if !haskey(plantDockDict, row[:Plant_code]*"__"*(ismissing(row[:Plant_dock]) ? "missing" : row[:Plant_dock]))
                nbPlantDocks = nbPlantDocks + 1
                plantDockDict[row[:Plant_code]*"__"*(ismissing(row[:Plant_dock]) ? "missing" : row[:Plant_dock])] = nbPlantDocks
            end
            
            
        end
    end
    

    
    TE_P = Matrix{Union{Float64, Missing}}(missing, nbPlannedTrucks, nbSuppliers)
    
    TL_P = Vector{Union{Float64, Missing}}(missing, nbPlannedTrucks)
    TW_P = Vector{Union{Float64, Missing}}(missing, nbPlannedTrucks)
    TH_P = Vector{Union{Float64, Missing}}(missing, nbPlannedTrucks)
    
    TKE_P = Matrix{Union{Float64, Missing}}(missing, nbPlannedTrucks, nbSupplierDocks)
    fill!(TKE_P, nbSupplierDocks)
    
    TGE_P = Matrix{Union{Float64, Missing}}(missing, nbPlannedTrucks, nbPlantDocks)
    fill!(TGE_P, nbPlantDocks)
    
    TDA_P = Vector{Union{Float64, Missing}}(missing, nbPlannedTrucks)
    
    TU_P = falses(nbPlannedTrucks, nbSuppliers)
    TP_P = falses(nbPlannedTrucks, nbPlants)
    TK_P = falses(nbPlannedTrucks, nbSupplierDocks)
    
    TG_P = falses(nbPlannedTrucks, nbPlantDocks)
    
    
    TR_P = falses(nbPlannedTrucks, nbItems) # TR is expanded, it will contain also only items which docks are stopped by by the truck
    
    
    # For each line of input_trucks:
    open(*(instancePath, "input_trucks.csv")) do input_trucksfile
        for row in CSV.File(input_trucksfile, normalizenames=true, delim=';', decimal=',', stripwhitespace=true, types=String)
            
            # Fill relevant truck information
            TE_P[truckDict[row[:Id_truck]], supplierDict[row[:Supplier_code]]] = parse(Float64, row[:Supplier_loading_order])
            TL_P[truckDict[row[:Id_truck]]] = parse(Float64, row[:Length])
            TW_P[truckDict[row[:Id_truck]]] = parse(Float64, row[:Width])
            TH_P[truckDict[row[:Id_truck]]] = parse(Float64, row[:Height])
            TKE_P[truckDict[row[:Id_truck]], supplierDockDict[row[:Supplier_code]*"__"*(ismissing(row[:Supplier_dock]) ? "missing" : row[:Supplier_dock])]] = parse(Float64, row[:Supplier_dock_loading_order])
            TGE_P[truckDict[row[:Id_truck]], plantDockDict[row[:Plant_code]*"__"*(ismissing(row[:Plant_dock]) ? "missing" : row[:Plant_dock])]] = parse(Float64, row[:Plant_dock_loading_order])
            TDA_P[truckDict[row[:Id_truck]]] = parse(Float64, row[:Arrival_time])
            TU_P[truckDict[row[:Id_truck]], supplierDict[row[:Supplier_code]]] = 1.0
            TP_P[truckDict[row[:Id_truck]], plantDict[row[:Plant_code]]] = 1.0


            TG_P[truckDict[row[:Id_truck]], plantDockDict[row[:Plant_code]*"__"*(ismissing(row[:Plant_dock]) ? "missing" : row[:Plant_dock])]] = 1.0
            TK_P[truckDict[row[:Id_truck]], supplierDockDict[row[:Supplier_code]*"__"*(ismissing(row[:Supplier_dock]) ? "missing" : row[:Supplier_dock])]] = 1.0
            
            
            # For each line of input trucks, retrieve Product code. For all items, get 
            # indices of item of same product code, and use it to fill TR
            product_code = row[:Product_code]
            for i in 1:nbItems
                TR_P[truckDict[row[:Id_truck]], i] = Iproduct_code[i] == product_code ? 1.0 : TR_P[truckDict[row[:Id_truck]], i]
            end
            # @debug begin
            #     if sum([Iproduct_code[i] == product_code ? 1.0 : 0.0 for i in 1:nbItems]) >= 1
            #         @debug product_code product_code
            #         @debug "[Iproduct_code[i] == product_code ? 1.0 : 0.0 for i in 1:nbItems]" [Iproduct_code[i] == product_code ? 1.0 : 0.0 for i in 1:10]
            #         @debug "TR_P[truckDict[row[:Id_truck]], :]" TR_P[truckDict[row[:Id_truck]], 1:10]
            #         sleep(20)
            #     end
            # end
        end
    end
    
    # IU = normalizeValues(input_itemsCSV[:Supplier_code])
    IU = falses(nbItems, nbSuppliers)
    # IP = normalizeValues(input_itemsCSV[:Plant_code])
    IP = falses(nbItems, nbPlants)
    # IK = normalizeValues(input_itemsCSV[:Supplier_dock])
    IK = falses(nbItems, nbSupplierDocks)
    # IPD = normalizeValues(input_itemsCSV[:Plant_dock])
    IPD = falses(nbItems, nbPlantDocks)
    # IS = normalizeValues(input_itemsCSV[:Stackability_code])
    IS = Vector{Union{Float64, Missing}}(missing, nbItems)
    
    IDL = Vector{Union{Float64, Missing}}(missing, nbItems)
    
    IDE = Vector{Union{Float64, Missing}}(missing, nbItems)
    
    itemDict = Dict{Int64, String}()

    stackabilitycodeDict = Dict{String, Float64}()
    nbstackabilitycodes = 0
    open(*(instancePath, "input_items.csv")) do input_itemsfile

        for (i, row) in enumerate(CSV.File(input_itemsfile, normalizenames=true, delim=';', decimal=',', stripwhitespace=true, types=String))
            IU[i, supplierDict[row[:Supplier_code]]] = 1.0
            IP[i, plantDict[row[:Plant_code]]] = 1.0
            IK[i, supplierDockDict[row[:Supplier_code]*"__"*(ismissing(row[:Supplier_dock]) ? "missing" : row[:Supplier_dock])]] = 1.0
            IPD[i, plantDockDict[row[:Plant_code]*"__"*(ismissing(row[:Plant_dock]) ? "missing" : row[:Plant_dock])]] = 1.0
            IDL[i] = parse(Float64, row[:Latest_arrival_time])
            IDE[i] = parse(Float64, row[:Earliest_arrival_time])
            if !haskey(stackabilitycodeDict, row[:Stackability_code])
                nbstackabilitycodes = nbstackabilitycodes + 1.0
                stackabilitycodeDict[row[:Stackability_code]] = nbstackabilitycodes
            end
            IS[i] = stackabilitycodeDict[row[:Stackability_code]]
            itemDict[i] = row[:Item_ident]
        end
    end
    # Expand TR with information about docks
    # For each truck, for each item, if the truck doesn't stop at the supplier & supplier dock of the item or 
    # it doesn't stop by the plant & plant dock of the item: replace with 0
    @debug sum(TR_P) sum(TR_P)

    # This loop might be useless or too restrictive, because it seems that the candidate list of 
    # each truck already includes this type of info 
    for t in 1:nbPlannedTrucks
        for i in 1:nbItems
            # if item `i` is compatible with truck `t` according to the candidate list
            if TR_P[t,i] == 1
                # if the supplier dock of i is not in the supplier docks of t
                # or the same but with plant docks
                # or the truck arrives after the latest arrival date allowed for i 
                if !*((IK[i,:] .<= TK_P[t,:])...) || !*((IPD[i,:] .<= TG_P[t,:])...) || TDA_P[t] > IDL[i]
                    # Remove i from compatible items for t
                    TR_P[t,i] = 0.0
                end
            end
        end
    end
    
    # @debug "" !*((IK[5,:] .<= TK_P[truckDict["P192711301"],:])...) || !*((IPD[5,:] .<= TG_P[truckDict["P192711301"],:])...) || TDA_P[truckDict["P192711301"]] > IDL[5]
    # @debug begin
    #     println(!*((IK[5,:] .<= TK_P[truckDict["P192711301"],:])...))
    #     println(!*((IPD[5,:] .<= TG_P[truckDict["P192711301"],:])...))
    #     println(TDA_P[truckDict["P192711301"]] > IDL[5])
    #     println("Item ident of 5", itemDict[5])
    # end
    # The total number of trucks could be nbPlannedTrucks * nbItems, but a 
    # smarter way would be to have
    # nbTrucks = sum(TR) # sum of number of candidate items for each truck
    # nbTrucks = nbPlannedTrucks * nbItems # Very bad idea (millions of trucks)
    nbTrucks = nbPlannedTrucks + sum([sum(TR_P[t, :]) > 0 ? sum(TR_P[t, :])-1 : 0 for t in 1:nbPlannedTrucks] )
    # @debug nbPlannedTrucks nbPlannedTrucks
    # @debug sum(TR_P) sum(TR_P)
    # @debug nbItems nbItems
    # @debug "" nbTrucks
    # @debug "nbPlannedTrucks * nbItems" nbPlannedTrucks * nbItems
    

    TE = Matrix{Union{Float64, Missing}}(missing, nbTrucks, nbSuppliers)
    
    TL = Vector{Union{Float64, Missing}}(missing, nbTrucks)
    TW = Vector{Union{Float64, Missing}}(missing, nbTrucks)
    TH = Vector{Union{Float64, Missing}}(missing, nbTrucks)
    
    TKE = Matrix{Union{Float64, Missing}}(missing, nbTrucks, nbSupplierDocks)
    fill!(TKE, nbSupplierDocks)
    
    TGE = Matrix{Union{Float64, Missing}}(missing, nbTrucks, nbPlantDocks)
    fill!(TGE, nbPlantDocks)
    
    TDA = Vector{Union{Float64, Missing}}(missing, nbTrucks)
    
    TU = falses(nbTrucks, nbSuppliers)
    TP = falses(nbTrucks, nbPlants)
    TK = falses(nbTrucks, nbSupplierDocks)
    
    TG = falses(nbTrucks, nbPlantDocks)
    
    
    TR = falses(nbTrucks, nbItems) # TR is expanded, it will contain also only items which docks are stopped by by the truck
    TID = Vector{Union{String, Missing}}(missing, nbTrucks)
    reverse_truckDict = Dict(value => key for (key, value) in truckDict)
    # For each planned truck
    let j = nbPlannedTrucks
        for p in 1:nbPlannedTrucks

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
            TID[p] = reverse_truckDict[p]

            TR[p, :] .= TR_P[p, :]
            
            # @debug p p
            # @debug "sum(TR_P[p,:])" sum(TR_P[p,:])
            # j = nbPlannedTrucks +sum(TR_P[1:p-1, :])-p
            # j = nbPlannedTrucks +sum(TR_P[1:p-1, :])-sum([sum(TR_P[k, :]) != 0 ? 1 : 0 for k in 1:p])

            # for each candidate items, add an extra truck
            for e in 1:sum(TR_P[p,:])-1
                j = j + 1
                # @debug begin
                #     if p in [3 8]
                #         println(e)
                #         println(j)
                #         sleep(1)
                #     end
                # end
                # Fill relevant truck information
                tail = "_E" * string(e)
                # if !haskey(truckDict, reverse_truckDict[p] * tail)
                #     truckDict[reverse_truckDict[p] * tail] = j
                #     reverse_truckDict[j] = reverse_truckDict[p] * tail
                # end
                TID[j] = reverse_truckDict[p] * tail
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



    open("tmp2.txt", "w") do io
        show(io, "text/plain", [string(reverse_truckDict[i], " : ", sum(TR_P[i,:])) for i in 1:size(TR_P)[1]])
    end

    
    open("tmp.txt", "w") do io
        show(io, "text/plain", TID)
    end

    # @info string("Nb of items with no truck available according to TR: ", sum(sum(TR_P[:, j]) == 0 ? 1 : 0 for j in 1:size(TR_P)[2]))
    if sum(sum(TR_P[:, j]) == 0 ? 1 : 0 for j in 1:size(TR_P)[2]) > 0
        @warn string("Nb of items with no truck available according to TR: ", sum(sum(TR_P[:, j]) == 0 ? 1 : 0 for j in 1:size(TR_P)[2]))
    end
    # @debug begin
    #     println([sum(TR_P[:, j]) == 0 ? string(j, "\n") : "" for j in 1:size(TR_P)[2]]...)
    # end

    nbStacks = nbItems

    model = Model(Cbc.Optimizer)
    @info "Creating variables..."
    @info "Adding zetaT..."
    @variable(model, zetaT[1:nbPlannedTrucks] >= 0)
    @info "Adding zetaE..."
    @variable(model, zetaE[1:nbTrucks - nbPlannedTrucks] >= 0)

    @info "Adding SS..."
    @variable(model, SS[1:nbStacks] >= 0)
    @info "Adding SP..."
    @variable(model, SP[1:nbStacks] >= 0)
    @info "Adding SK..."
    @variable(model, SK[1:nbStacks] >= 0)
    @info "Adding SPD..."
    @variable(model, SPD[1:nbStacks] >= 0)
    @info "Adding SU..."
    @variable(model, SU[1:nbStacks] >= 0)
    @info "Adding SO..."
    @variable(model, SO[1:nbStacks] >= 0)
    @info "Adding SXe..."
    @variable(model, SXe[1:nbStacks] >= 0)
    @info "Adding SXo..."
    @variable(model, SXo[1:nbStacks] >= 0)
    @info "Adding SYe..."
    @variable(model, SYe[1:nbStacks] >= 0)
    @info "Adding SYo..."
    @variable(model, SYo[1:nbStacks] >= 0)
    @info "Adding SZe..."
    @variable(model, SZe[1:nbStacks] >= 0)
    @info "Adding betaM..."
    @variable(model, betaM[1:nbStacks] >= 0)
    @info "Adding betaP..."
    @variable(model, betaP[1:nbStacks] >= 0)
    @info "Adding nu..."
    @variable(model, nu[1:nbStacks] >= 0)
    @info "Adding tau..."
    @variable(model, tau[1:nbStacks] >= 0)
    @info "Adding phi..."
    @variable(model, phi[1:nbStacks] >= 0)
    @info "Adding SG..."
    @variable(model, SG[1:nbStacks, 1:nbPlants] >= 0)

    @info "Adding TI..."
    @variable(model, TI[1:nbTrucks, 1:nbItems], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding R..."
    @variable(model, R[1:nbItems, 1:nbSuppliers], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding Theta..."
    @variable(model, Theta[1:nbItems, 1:nbSuppliers], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding S..."
    @variable(model, S[1:nbStacks, 1:nbItems], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding Z..."
    @variable(model, Z[1:nbStacks, 1:nbItems], lower_bound = 0, upper_bound = 1, Bin)

    # Omega is too big to handle, even in small instances TODO
    @debug nbStacks nbStacks
    @debug nbTrucks nbTrucks
    @debug nbItems nbItems
    @debug "nbStacks * nbTrucks * nbItems" nbStacks * nbTrucks * nbItems
    @info "Adding Omega..."
    # @variable(model, Omega[1:nbStacks, 1:nbTrucks, 1:nbItems], lower_bound = 0, upper_bound = 1, container=Array, Bin)

    # TODO this is taking way too long
    for t in 1:nbTrucks
        @info string("Adding Omega[", t, "]...")
        model[Symbol("Omega[", t, "]")] = @variable(model, [1:nbStacks, 1:nbItems], lower_bound = 0, upper_bound = 1, container=Array, Bin)
    end

    @info "Adding ST..."
    @variable(model, ST[1:nbStacks, 1:nbTrucks], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding IOV..."
    @variable(model, IOV[1:nbItems], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding mu..."
    @variable(model, mu[1:nbStacks], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding eta..."
    @variable(model, eta[1:nbStacks], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding xi..."
    @variable(model, xi[1:nbStacks], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding chi..."
    @variable(model, chi[1:nbStacks], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding r..."
    @variable(model, r[1:nbStacks], lower_bound = 0, upper_bound = 1, Bin)

    @info "Adding sigma1..."
    @variable(model, sigma1[1:nbStacks], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding sigma2..."
    @variable(model, sigma2[1:nbStacks], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding sigma3..."
    @variable(model, sigma3[1:nbStacks], lower_bound = 0, upper_bound = 1, Bin)

    @info "Adding Psi..."
    @variable(model, Psi[1:nbStacks, 1:nbTrucks], lower_bound = 0, Int)
    @info "Adding Q..."
    @variable(model, Q[1:nbStacks, 1:nbItems], lower_bound = 0, Int)
    @info "Adding H..."
    @variable(model, H[1:nbStacks, 1:nbItems], lower_bound = 0, Int)
    @info "Adding V..."
    @variable(model, V[1:nbStacks, 1:nbItems], lower_bound = 0, Int)
    @info "Adding W..."
    @variable(model, W[1:nbStacks, 1:nbItems], lower_bound = 0, Int)
    @info "Adding Gl..."
    @variable(model, Gl[1:nbStacks, 1:nbItems], lower_bound = 0, Int)
    @info "Adding Gr..."
    @variable(model, Gr[1:nbStacks, 1:nbItems], lower_bound = 0, Int)

    @info "Adding lambda..."
    @variable(model, lambda[[1:nbStacks * (nbStacks+1)/2]], lower_bound = 0, Bin)
    
    @info "Computing parameters..."
    # MI4 = [[max(TU[:, j]...) for j in 1:first(size(TU[1, :]))] for j in 1:first(size(TI[:, 1]))]
    MI4 = map(Int, (max(TU[:, i]...) for i = 1:first(size(TU[1, :])), j = 1:first(size(TI[:, 1]))))

    MZ = max(IS...) + 1.0

    MQ = max(IU...) + 1.0

    MH = max(IP...) + 1.0

    MV = max(IK...) + 1.0

    MW = max(IPD...) + 1.0

    Gr = max(SO...) + 1.0

    Gl = max(IO...) + 1.0

    Meta = nbTrucks * Mtau/10

    MTE = max(TE...)

    MTKE = max(TKE...)

    MTGE = max(TGE...)

    MTW = max(TW...) + 1.0

    MPsi = nbItems

    Mlambda = 2.0 * TL .+ 1.0 # TODO check if Mlambda big enough

    SZo = 0.0

    MST = 2.0

    MOmega = 2.0

    Mmu = 2.0

    Mtau = 10.0

    Xi1 = falses(convert(Int, nbStacks*(nbStacks+1)/2), nbStacks)
    Xi2 = falses(convert(Int, nbStacks*(nbStacks+1)/2), nbStacks)

    epsilon = 0.001


    fillXi1!(Xi1)
    fillXi2!(Xi2)

    @info "Adding constraints..."
    @info "Adding cZetaT1..."
    @constraint(model, cZetaT1, -zetaT >= -ones(size(zetaT)[1]))
    
    @info "Adding cZetaT2..."
    @constraint(model, cZetaT2, -zetaT >= -TIT * ones(size(zetaT)[1]))

    @info "Adding cZetaE1..."
    @constraint(model, cZetaE1, -zetaT >= -ones(size(zetaT)[1]))
    @info "Adding cZetaE2..."
    @constraint(model, cZetaE2, -zetaT >= -TIE * ones(size(zetaT)[1]))

    @info "Adding cTI_TR..."
    @constraint(model, cTI_TR, TI <= TR)

    @info "Adding cTI_1_1..."
    @constraint(model, cTI_1_1, transpose(TI) * ones(size(transpose(TI))[1]) <= ones(size(transpose(TI))[1]))

    @info "Adding cTI_TP_IP..."
    @constraint(model, cTI_TP_IP, transpose(TI) * TP == IP)

    @info "Adding cR_Theta_MI4..."
    @constraint(model, cR_Theta_MI4, R <= Theta * MI4)

    @info "Adding cR_TI_TU..."
    @constraint(model, cR_TI_TU, -MI4*(1-Theta) <= R - (transpose(TI) * TU) <= MI4* (1-Theta))

    @info "Adding cR_1_IU..."
    @constraint(model, cR_1_IU, R * ones(size(R)[1]) == IU)

    @info "Adding cTheta_1_1..."
    @constraint(model, cTheta_1_1, Theta * ones(size(Theta)[1]) >= ones(size(Theta)[1]))

    @info "Adding cIDE_TI_TDA_IDL..."
    @constraint(model, cIDE_TI_TDA_IDL, IDE <= transpose(TI) * TDA <= IDL)

    @info "Adding cZ_S_MZ..."
    @constraint(model, cZ_S_MZ, Z <= S * MZ)

    @info "Adding cPsi_Omega..."
    @constraint(model, cPsi_Omega, Psi == hcat([Omega[:,i,:]*ones(size(Omega)[1]) for i in 1:size(Omega)[2]]...))

    for j in 1:size(Omega)[2]
        @info "Adding cOmega_S_MOmega..."
        # @constraint(model, cOmega_S_MOmega, Omega[:, j, :] <= S * MOmega)
        @constraint(model, cOmega_S_MOmega, model(Symbol("Omega[", j, "]")) <= S * MOmega)
    end
    # for i in 1:size(Omega)[1]
    #     @info "Adding cOmega_TI_S..."
    #     @constraint(model, cOmega_TI_S, -(1-S)*MOmega <= Omega[i,:,:] - transpose(TI) <= (1 - S)*MOmega)
    # end

    @info "Adding cPsi_ST_MPsi..."
    @constraint(model, cPsi_ST_MPsi, Psi <= St * MPsi)

    @info "Adding cPsi_S_ST..."
    @constraint(model, cPsi_S_ST, -(1-ST)*MPsi <= Psi - hcat([S * ones(size(S)[1]) for i in 1:size(Psi)[2]]) <= (1-ST)*MPsi)

    @info "Adding cS_IS_Z..."
    @constraint(model, cS_IS_Z, S * IS == Z * ones(size(Z)[1]))

    @info "Adding cZ_SS_S..."
    @constraint(model, cZ_SS_S, -MZ * (1-S) <= Z - hcat([SS for i in 1:size(Z)[2]]) <= MZ * (1-S))

    @info "Adding cQ_S..."
    @constraint(model, cQ_S, Q <= S * MQ)

    @info "Adding cS_IU_Q..."
    @constraint(model, cS_IU_Q, S * IU == Q * ones(size(Q)[1]))

    @info "Adding cQ_SU_S..."
    @constraint(model, cQ_SU_S, -MQ * (1-S) <= Q - hcat([SU for i in 1:size(Q)[2]]) <= MQ * (1-S))


    @info "Adding cH_S..."
    @constraint(model, cH_S, H <= S * MH)

    @info "Adding cS_IP_H..."
    @constraint(model, cS_IP_H, S * IP == H * ones(size(H)[1]))

    @info "Adding cQ_SP_S..."
    @constraint(model, cQ_SP_S, -MH * (1-S) <= H - hcat([SP for i in 1:size(H)[2]]) <= MH * (1-S))


    @info "Adding cV_S..."
    @constraint(model, cV_S, V <= S * MV)

    @info "Adding cS_IK_V..."
    @constraint(model, cS_IK_V, S * IK == V * ones(size(V)[1]))

    @info "Adding cV_SK_S..."
    @constraint(model, cV_SK_S, -MV * (1-S) <= V - hcat([SK for i in 1:size(V)[2]]) <= MV * (1-S))


    @info "Adding cW_S..."
    @constraint(model, cW_S, W <= S * MW)

    @info "Adding cS_IPD_W..."
    @constraint(model, cS_IPD_W, S * IPD == W * ones(size(W)[1]))

    @info "Adding cW_SPD_S..."
    @constraint(model, cW_SPD_S, -MW * (1-S) <= W - hcat([SPD for i in 1:size(W)[2]]) <= MW * (1-S))


    @info "Adding cGl_S..."
    @constraint(model, cGl_S, Gl <= S * MG)
    @info "Adding cGr_S..."
    @constraint(model, cGr_S, Gr <= S * MG)

    @info "Adding Gl..."
    @constraint(model, Gl * ones(size(Gl)[1]) == Gr * ones(size(Gr)[1]))

    @info "Adding cGr_SO_S..."
    @constraint(model, cGr_SO_S, -MG * (1-S) <= Gr - hcat([SO for i in 1:size(Gr)[2]]) <= MG * (1-S))
    @info "Adding cGl_IOV_S..."
    @constraint(model, cGl_IOV_S, -MG * (1-S) <= Gl - hcat([IOV for i in 1:size(Gl)[2]]) <= MG * (1-S))

    @info "Adding cSXe_SXo_SL_SO..."
    @constraint(model, cSXe_SXo_SL_SO, SXe - SXo == SL + SO * MTL)
    @info "Adding cSYe_SYo_SW_SO..."
    @constraint(model, cSYe_SYo_SW_SO, SYe - SYo == SW + SO * MTW)

    @info "Adding cSXe_SXo_SW_SO..."
    @constraint(model, cSXe_SXo_SW_SO, SXe - SXo == SW + (1 - SO) * MTW)
    @info "Adding cSYe_SYo_SL_SO..."
    @constraint(model, cSYe_SYo_SL_SO, SYe - SYo == SL + (1 - SO) * MTL)

    @info "Adding cSZe_S_IH..."
    @constraint(model, cSZe_S_IH, SZe == S* IH)

    @info "Adding cSXe_ST_TL..."
    @constraint(model, cSXe_ST_TL, SXe <= ST * TL)
    @info "Adding cSYe_ST_TW..."
    @constraint(model, cSYe_ST_TW, SYe <= ST * TW)
    @info "Adding cSZe_ST_TH..."
    @constraint(model, cSZe_ST_TH, SZe <= ST * TH)

    @info "Adding cSXo_SXo..."
    @constraint(model, cSXo_SXo, vcat(hcat([1], falses(1, size(SXo)[1]-1)), I(size(SXo)[1])) * SXo <= SXo)

    @info "Adding cXi2SXo_Xi1SXe_betaM_betaP..."
    @constraint(model, cXi2SXo_Xi1SXe_betaM_betaP, Xi2 * SXo - Xi1 * SXe - betaM + betaP == -epsilon)

    @info "Adding cbetaM_lambda..."
    @constraint(model, cbetaM_lambda, betaM <= lambda * Mlambda)

    @info "Adding betaP..."
    @constraint(model, betaP <= (1-lambda)*Mlambda)

    @info "Adding cmu_betaM..."
    @constraint(model, cmu_betaM, (1-mu) <= betaM * Mmu)

    @info "Adding cXi2ST_Xi1ST_nu..."
    @constraint(model, cXi2ST_Xi1ST_nu, (Xi2 * ST - Xi1 * ST) * diagm([i for i in 1:size(nu)[2]]) == nu)

    @info "Adding ctau_phi_nu..."
    @constraint(model, ctau_phi_nu, tau - phi <= (nu - nbTrucks) * Mtau)

    @info "Adding ctau_nu..."
    @constraint(model, ctau_nu, tau >= (nu-nbTrucks)*Mtau/10)

    @info "Adding ctau_eta..."
    @constraint(model, ctau_eta, tau <= eta * Meta)

    @info "Adding cphi_eta..."
    @constraint(model, cphi_eta, phi <= (1 - eta)*Meta)

    @info "Adding cXi1SYe_Xi2SYo..."
    @constraint(model, cXi1SYe_Xi2SYo, Xi1 * SYe <= Xi2 * SYo + xi * MTW + (tau + phi) * MTW + (1 - mu) * MTW)
    @info "Adding cXi2SYe_Xi1SYo..."
    @constraint(model, cXi2SYe_Xi1SYo, Xi2 * SYe <= Xi1 * SYo + (1-xi) * MTW + (tau + phi) * MTW + (1 - mu) * MTW)

    @info "Adding cXi1SU_Xi2SU..."
    @constraint(model, cXi1SU_Xi2SU, Xi1 * SU * TE <= Xi2 * SU * TE + (tau + phi) * MTE)

    @info "Adding cXi1SU_Xi2SU_chi..."
    @constraint(model, cXi1SU_Xi2SU_chi, Xi1SU - Xi2SU >= chi * epsilon - r*MTE - (tau + phi) * MTE - (1 - sigma1) * MTE)

    @info "Adding cXi2SU_Xi1SU_chi..."
    @constraint(model, cXi2SU_Xi1SU_chi, Xi2SU - Xi1SU >= (1-chi) * epsilon - r*MTE - (tau + phi) * MTE - (1 - sigma1) * MTE)

    @info "Adding cXi2SK_Xi1SK..."
    @constraint(model, cXi2SK_Xi1SK, Xi2*SK*TKE >= Xi1*SK*TKE - (1 - r) * MTKE - (tau + phi) * MTKE)

    @info "Adding cXi1SK_Xi2SK_chi..."
    @constraint(model, cXi1SK_Xi2SK_chi, Xi1*SK*TKE - Xi2*SK*TKE >= chi*epsilon - (tau + phi)*MTKE - (1 - sigma2)*MTKE)
    @info "Adding cXi2SK_Xi1SK_chi..."
    @constraint(model, cXi2SK_Xi1SK_chi, Xi2*SK*TKE - Xi1*SK*TKE >= (1-chi)*epsilon - (tau + phi)*MTKE - (1 - sigma2)*MTKE)

    @info "Adding cXi2SG_Xi1SG..."
    @constraint(model, cXi2SG_Xi1SG, Xi2*SG*TGE >= Xi1*SG*TGE - (tau + phi)*MTGE - (1 - sigma3) * MTGE)

    @info "Adding csigma1_sigma2_sigma3..."
    @constraint(model, csigma1_sigma2_sigma3, sigma1 + sigma2 + sigma3 >= 1)


    @info "Displaying model..."
    display(model)
end
    @allocated main()

end
