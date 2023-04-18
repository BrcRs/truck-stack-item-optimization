module TSIInstanceLoader

using LinearAlgebra
using JuMP
using Cbc
using CSV
using FilePaths
using Logging

logger = ConsoleLogger(stdout, Logging.Debug)

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

function main()
    instancePath = "Instances/AS/"

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

    # input_itemsCSV = CSV.File(input_itemsfile, normalizenames=true, delim=';', decimal=',', stripwhitespace=true)
    # input_parametersCSV = CSV.File(input_parametersfile, normalizenames=true, delim=';', decimal=',', stripwhitespace=true)
    # input_trucksCSV = CSV.File(input_trucksfile, normalizenames=true, delim=';', decimal=',', stripwhitespace=true)
    # input_parametersCSV = CSV.read(open(Path(instancePath, "input_parameters.csv")), normalizenames=true, delim=";", decimal=",", stripwhitespace=true)
    # input_trucksCSV = CSV.read(open(Path(instancePath, "input_trucks.csv")), normalizenames=true, delim=";", decimal=",", stripwhitespace=true)
    nbItems = 0
    Iproduct_code = Vector{String}()
    open(*(instancePath, "input_items.csv")) do input_itemsfile
        for row in CSV.Rows(input_itemsfile, normalizenames=true, delim=';', decimal=',', stripwhitespace=true)
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
        for row in CSV.Rows(input_trucksfile, normalizenames=true, delim=';', decimal=',', stripwhitespace=true)
            if !haskey(truckDict, row[:Id_truck])
                nbPlannedTrucks = nbPlannedTrucks + 1
                truckDict[row[:Id_truck]] = nbPlannedTrucks
            end
            if !haskey(supplierDict, row[:Supplier_code])
                nbSuppliers = nbSuppliers + 1
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
        for row in CSV.Rows(input_trucksfile, normalizenames=true, delim=';', decimal=',', stripwhitespace=true)
            
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
            TR_P[truckDict[row[:Id_truck]], :] .= [Iproduct_code[i] == product_code ? 1.0 : 0.0 for i in 1:nbItems]
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
    
    stackabilitycodeDict = Dict{String, Float64}()
    nbstackabilitycodes = 0
    open(*(instancePath, "input_items.csv")) do input_itemsfile

        for (i, row) in enumerate(CSV.Rows(input_itemsfile, normalizenames=true, delim=';', decimal=',', stripwhitespace=true))
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
        end
    end
    # Expand TR with information about docks
    # For each truck, for each item, if the truck doesn't stop at the supplier & supplier dock of the item or 
    # it doesn't stop by the plant & plant dock of the item: replace with 0
    @debug sum(TR_P) sum(TR_P)

    for t in 1:nbPlannedTrucks
        for i in 1:nbItems
            if TR_P[t,i] == 1
                if !*((IK[i,:] .<= TK_P[t,:])...) || !*((IPD[i,:] .<= TG_P[t,:])...) || TDA_P[t] > IDL[i]
                    TR_P[t,i] = 0.0
                end
            end
        end
    end
    
    
    # The total number of trucks could be nbPlannedTrucks * nbItems, but a 
    # smarter way would be to have
    # nbTrucks = sum(TR) # sum of number of candidate items for each truck
    # nbTrucks = nbPlannedTrucks * nbItems # Very bad idea (millions of trucks)
    nbTrucks = sum(TR_P)
    @debug nbPlannedTrucks nbPlannedTrucks
    @debug sum(TR_P) sum(TR_P)
    @debug nbItems nbItems
    @debug "" nbTrucks
    @debug "nbPlannedTrucks * nbItems" nbPlannedTrucks * nbItems
    

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
    

    # For each planned truck
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
        
        TR[p, :] .= TR_P[p, :]
        
        # for each candidate items, add an extra truck
        for e in 1:sum(TR_P[p,:])
            j = nbPlannedTrucks +e +sum(TR_P[1:p-1, :])
            # Fill relevant truck information
            tail = "_E" * convert(String, e)
            if !haskey(truckDict, row[:Id_truck] * tail)
                truckDict[row[:Id_truck] * tail] = truckDict[row[:Id_truck]] + e
            end
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
    end
    
        
    model = Model(Cbc.Optimizer)

    @variable(model, zetaT[1:nbPlannedTrucks] >= 0)
    @variable(model, zetaE[1:nbTrucks - nbPlannedTrucks] >= 0)

    @variable(model, SS[1:nbStacks] >= 0)
    @variable(model, SP[1:nbStacks] >= 0)
    @variable(model, SK[1:nbStacks] >= 0)
    @variable(model, SPD[1:nbStacks] >= 0)
    @variable(model, SU[1:nbStacks] >= 0)
    @variable(model, SO[1:nbStacks] >= 0)
    @variable(model, SXe[1:nbStacks] >= 0)
    @variable(model, SXo[1:nbStacks] >= 0)
    @variable(model, SYe[1:nbStacks] >= 0)
    @variable(model, SYo[1:nbStacks] >= 0)
    @variable(model, SZe[1:nbStacks] >= 0)
    @variable(model, betaM[1:nbStacks] >= 0)
    @variable(model, betaP[1:nbStacks] >= 0)
    @variable(model, nu[1:nbStacks] >= 0)
    @variable(model, tau[1:nbStacks] >= 0)
    @variable(model, phi[1:nbStacks] >= 0)
    @variable(model, SG[1:nbStacks, 1:nbPlants] >= 0)

    @variable(model, TI[1:nbTrucks, 1:sum(TR_P)], lower_bound = 0, upper_bound = 1)
    @variable(model, R[1:nbItems, 1:nbSuppliers], lower_bound = 0, upper_bound = 1)
    @variable(model, Theta[1:nbItems, 1:nbSuppliers], lower_bound = 0, upper_bound = 1)
    @variable(model, S[1:nbStacks, 1:nbItems], lower_bound = 0, upper_bound = 1)
    @variable(model, Z[1:nbStacks, 1:nbItems], lower_bound = 0, upper_bound = 1)
    @variable(model, Omega[1:nbStacks, 1:nbTrucks, 1:nbItems], lower_bound = 0, upper_bound = 1)
    @variable(model, ST[1:nbStacks, 1:nbTrucks], lower_bound = 0, upper_bound = 1)
    @variable(model, IOV[1:nbItems], lower_bound = 0, upper_bound = 1)
    @variable(model, mu[1:nbStacks], lower_bound = 0, upper_bound = 1)
    @variable(model, eta[1:nbStacks], lower_bound = 0, upper_bound = 1)
    @variable(model, xi[1:nbStacks], lower_bound = 0, upper_bound = 1)
    @variable(model, chi[1:nbStacks], lower_bound = 0, upper_bound = 1)
    @variable(model, r[1:nbStacks], lower_bound = 0, upper_bound = 1)

    @variable(model, sigma1[1:nbStacks], lower_bound = 0, upper_bound = 1)
    @variable(model, sigma2[1:nbStacks], lower_bound = 0, upper_bound = 1)
    @variable(model, sigma3[1:nbStacks], lower_bound = 0, upper_bound = 1)

    @variable(model, Psi[1:nbStacks, 1:nbTrucks], lower_bound = 0)
    @variable(model, Q[1:nbStacks, 1:nbItems], lower_bound = 0)
    @variable(model, H[1:nbStacks, 1:nbItems], lower_bound = 0)
    @variable(model, V[1:nbStacks, 1:nbItems], lower_bound = 0)
    @variable(model, W[1:nbStacks, 1:nbItems], lower_bound = 0)
    @variable(model, Gl[1:nbStacks, 1:nbItems], lower_bound = 0)
    @variable(model, Gr[1:nbStacks, 1:nbItems], lower_bound = 0)

    @variable(model, lambda[[1:nbStacks * (nbStacks+1)/2]], lower_bound = 0)

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

    @constraint(model, cZetaT1, -zetaT >= -ones(size(zetaT)[1]))
    @constraint(model, cZetaT2, -zetaT >= -TIT * ones(size(zetaT)[1]))

    @constraint(model, cZetaE1, -zetaT >= -ones(size(zetaT)[1]))
    @constraint(model, cZetaE2, -zetaT >= -TIE * ones(size(zetaT)[1]))

    @constraint(model, cTI_TR, TI <= TR)

    @constraint(model, cTI_1_1, transpose(TI) * ones(size(transpose(TI))[1]) <= ones(size(transpose(TI))[1]))

    @constraint(model, cTI_TP_IP, transpose(TI) * TP == IP)

    @constraint(model, cR_Theta_MI4, R <= Theta * MI4)

    @constraint(model, cR_TI_TU, -MI4*(1-Theta) <= R - (transpose(TI) * TU) <= MI4* (1-Theta))

    @constraint(model, cR_1_IU, R * ones(size(R)[1]) == IU)

    @constraint(model, cTheta_1_1, Theta * ones(size(Theta)[1]) >= ones(size(Theta)[1]))

    @constraint(model, cIDE_TI_TDA_IDL, IDE <= transpose(TI) * TDA <= IDL)

    @constraint(model, cZ_S_MZ, Z <= S * MZ)

    @constraint(model, cPsi_Omega, Psi == hcat([Omega[:,i,:]*ones(size(Omega)[1]) for i in 1:size(Omega)[2]]...))

    for j in 1:size(Omega)[2]
        @constraint(model, cOmega_S_MOmega, Omega[:, j, :] <= S * MOmega)
    end
    for i in 1:size(Omega)[1]
        @constraint(model, cOmega_TI_S, -(1-S)*MOmega <= Omega[i,:,:] - transpose(TI) <= (1 - S)*MOmega)
    end

    @constraint(model, cPsi_ST_MPsi, Psi <= St * MPsi)

    @constraint(model, cPsi_S_ST, -(1-ST)*MPsi <= Psi - hcat([S * ones(size(S)[1]) for i in 1:size(Psi)[2]]) <= (1-ST)*MPsi)

    @constraint(model, cS_IS_Z, S * IS == Z * ones(size(Z)[1]))

    @constraint(model, cZ_SS_S, -MZ * (1-S) <= Z - hcat([SS for i in 1:size(Z)[2]]) <= MZ * (1-S))

    @constraint(model, cQ_S, Q <= S * MQ)

    @constraint(model, cS_IU_Q, S * IU == Q * ones(size(Q)[1]))

    @constraint(model, cQ_SU_S, -MQ * (1-S) <= Q - hcat([SU for i in 1:size(Q)[2]]) <= MQ * (1-S))


    @constraint(model, cH_S, H <= S * MH)

    @constraint(model, cS_IP_H, S * IP == H * ones(size(H)[1]))

    @constraint(model, cQ_SP_S, -MH * (1-S) <= H - hcat([SP for i in 1:size(H)[2]]) <= MH * (1-S))


    @constraint(model, cV_S, V <= S * MV)

    @constraint(model, cS_IK_V, S * IK == V * ones(size(V)[1]))

    @constraint(model, cV_SK_S, -MV * (1-S) <= V - hcat([SK for i in 1:size(V)[2]]) <= MV * (1-S))


    @constraint(model, cW_S, W <= S * MW)

    @constraint(model, cS_IPD_W, S * IPD == W * ones(size(W)[1]))

    @constraint(model, cW_SPD_S, -MW * (1-S) <= W - hcat([SPD for i in 1:size(W)[2]]) <= MW * (1-S))


    @constraint(model, cGl_S, Gl <= S * MG)
    @constraint(model, cGr_S, Gr <= S * MG)

    @constraint(model, Gl * ones(size(Gl)[1]) == Gr * ones(size(Gr)[1]))

    @constraint(model, cGr_SO_S, -MG * (1-S) <= Gr - hcat([SO for i in 1:size(Gr)[2]]) <= MG * (1-S))
    @constraint(model, cGl_IOV_S, -MG * (1-S) <= Gl - hcat([IOV for i in 1:size(Gl)[2]]) <= MG * (1-S))

    @constraint(model, cSXe_SXo_SL_SO, SXe - SXo == SL + SO * MTL)
    @constraint(model, cSYe_SYo_SW_SO, SYe - SYo == SW + SO * MTW)

    @constraint(model, cSXe_SXo_SW_SO, SXe - SXo == SW + (1 - SO) * MTW)
    @constraint(model, cSYe_SYo_SL_SO, SYe - SYo == SL + (1 - SO) * MTL)

    @constraint(model, cSZe_S_IH, SZe == S* IH)

    @constraint(model, cSXe_ST_TL, SXe <= ST * TL)
    @constraint(model, cSYe_ST_TW, SYe <= ST * TW)
    @constraint(model, cSZe_ST_TH, SZe <= ST * TH)

    @constraint(model, cSXo_SXo, vcat(hcat([1], falses(1, size(SXo)[1]-1)), I(size(SXo)[1])) * SXo <= SXo)

    @constraint(model, cXi2SXo_Xi1SXe_betaM_betaP, Xi2 * SXo - Xi1 * SXe - betaM + betaP == -epsilon)

    @constraint(model, cbetaM_lambda, betaM <= lambda * Mlambda)

    @constraint(model, betaP <= (1-lambda)*Mlambda)

    @constraint(model, cmu_betaM, (1-mu) <= betaM * Mmu)

    @constraint(model, cXi2ST_Xi1ST_nu, (Xi2 * ST - Xi1 * ST) * diagm([i for i in 1:size(nu)[2]]) == nu)

    @constraint(model, ctau_phi_nu, tau - phi <= (nu - nbTrucks) * Mtau)

    @constraint(model, ctau_nu, tau >= (nu-nbTrucks)*Mtau/10)

    @constraint(model, ctau_eta, tau <= eta * Meta)

    @constraint(model, cphi_eta, phi <= (1 - eta)*Meta)

    @constraint(model, cXi1SYe_Xi2SYo, Xi1 * SYe <= Xi2 * SYo + xi * MTW + (tau + phi) * MTW + (1 - mu) * MTW)
    @constraint(model, cXi2SYe_Xi1SYo, Xi2 * SYe <= Xi1 * SYo + (1-xi) * MTW + (tau + phi) * MTW + (1 - mu) * MTW)

    @constraint(model, cXi1SU_Xi2SU, Xi1 * SU * TE <= Xi2 * SU * TE + (tau + phi) * MTE)

    @constraint(model, cXi1SU_Xi2SU_chi, Xi1SU - Xi2SU >= chi * epsilon - r*MTE - (tau + phi) * MTE - (1 - sigma1) * MTE)

    @constraint(model, cXi2SU_Xi1SU_chi, Xi2SU - Xi1SU >= (1-chi) * epsilon - r*MTE - (tau + phi) * MTE - (1 - sigma1) * MTE)

    @constraint(model, cXi2SK_Xi1SK, Xi2*SK*TKE >= Xi1*SK*TKE - (1 - r) * MTKE - (tau + phi) * MTKE)

    @constraint(model, cXi1SK_Xi2SK_chi, Xi1*SK*TKE - Xi2*SK*TKE >= chi*epsilon - (tau + phi)*MTKE - (1 - sigma2)*MTKE)
    @constraint(model, cXi2SK_Xi1SK_chi, Xi2*SK*TKE - Xi1*SK*TKE >= (1-chi)*epsilon - (tau + phi)*MTKE - (1 - sigma2)*MTKE)

    @constraint(model, cXi2SG_Xi1SG, Xi2*SG*TGE >= Xi1*SG*TGE - (tau + phi)*MTGE - (1 - sigma3) * MTGE)

    @constraint(model, csigma1_sigma2_sigma3, sigma1 + sigma2 + sigma3 >= 1)



    display(model)

end
    main()

end
