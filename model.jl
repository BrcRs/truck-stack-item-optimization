using JuMP

include("instance_loader.jl")
include("matrix_ops.jl")

mutable struct Subproblem
    t::Integer
    problem::TSIProblem
    submodel::Model
    optimizer
    nbstacks::Integer
end

mutable struct TSIProblem
    model::Model
    obj_dict::Dict{Symbol,Any}

end

object_dictionary(pb::TSIProblem) = pb.obj_dict

function unregister(pb::TSIProblem, key::Symbol)
    return JuMP.unregister(pb.model, key)
end

function Base.getindex(p::TSIProblem, name::Symbol)
    obj_dict = object_dictionary(p)
    if !haskey(obj_dict, name)
        throw(KeyError(name))
    end
    return obj_dict[name]
end

function Base.setindex!(problem::TSIProblem, value, name::Symbol)
    return object_dictionary(problem)[name] = value
end

function unregister(problem::TSIProblem, key::Symbol)
    delete!(object_dictionary(problem), key)
    return
end
function Base.haskey(problem::TSIProblem, name::Symbol)
    return haskey(object_dictionary(problem), name)
end


function solve!(problem::TSIProblem)
    # 1. Make first solution by distributing items in planned trucks
    # Allocating TIbar
    # TIbar = vcat(problem[:TR_P], falses(problem[:nbtrucks] - problem[:nbplannedtrucks], problem[:nbitems]))

    TIbar = falses(nbtrucks, nbitems)
    # For each item, assign a random allowed truck
    for i in 1:nbitems
        icandidates = findall((x) -> x == 1, TIbar[:, i])
        TIbar[rand(icandidates)] = 1
    end

    # 2. Instanciate as many subproblems than there are trucks. 
    # In each subproblem, only create constraints related to the corresponding truck.
    # Make predetermined stacks
    subproblems = [Subproblem(t, problem) for t in 1:nbtrucks]
    # while TI[t, k+1] - sum([pi[t] * TI[t, k+1] for t in bold_T]) != 0
    #     for t in bold_T
    #         # Solve the deterministic minimization problem for truck t with
    #         # a penalization + kappa[t, k] * (TI[t, k+1] - TIbar[k])
    #         # and obtain optimal first decision TI[t, k+1]
    #         ...
    #     end
    #     # Update the mean first decisions:
    #     TIbar[k+1] = sum([pi[t] * TI[t, k+1] for t in bold_T])

    #     # Update the multipliers by
    #     for t in bold_T
    #         kappa[t, k+1] = kappa[t, k] + delta * (TI[t, k+1] - TIbar[k+1])
    #     end
end

function TSIProblem(optimizer, instancepath)
    item_productcodes, truckdict, supplierdict, supplierdockdict, plantdict, 
    plantdockdict, nbplannedtrucks, nbitems, nbsuppliers, nbsupplierdocks, nbplants, 
    nbplantdocks, TE_P, TL_P, TW_P, TH_P, TKE_P, TGE_P, TDA_P, TU_P, TP_P, TK_P, 
    TG_P, TR_P, IU, IP, IK, IPD, IS, IDL, IDE, stackabilitycodedict, nbtrucks, TE, 
    TL, TW, TH, TKE, TGE, TDA, TU, TP, TK, TG, TR, TID, reverse_truckdict,
    costinventory, costtransportation, costextratruck, timelimit = loadinstance(instancepath)
    
    nbstacks = nbitems

    return TSIProblem(optimizer, 
    item_productcodes, truckdict, supplierdict, supplierdockdict, plantdict, 
    plantdockdict, nbplannedtrucks, nbitems, nbsuppliers, nbsupplierdocks, nbplants, 
    nbplantdocks, nbstacks, TE_P, TL_P, TW_P, TH_P, TKE_P, TGE_P, TDA_P, TU_P, TP_P, TK_P, 
    TG_P, TR_P, IU, IP, IK, IPD, IS, IDL, IDE, stackabilitycodedict, nbtrucks, TE, 
    TL, TW, TH, TKE, TGE, TDA, TU, TP, TK, TG, TR, TID, reverse_truckdict,  
    costinventory, costtransportation, costextratruck, timelimit
    )

end

getmodel(pb::TSIProblem) = pb.model

function TSIProblem(optimizer, 
    item_productcodes, truckdict, supplierdict, supplierdockdict, plantdict, 
    plantdockdict, nbplannedtrucks, nbitems, nbsuppliers, nbsupplierdocks, nbplants, 
    nbplantdocks, nbstacks, TE_P, TL_P, TW_P, TH_P, TKE_P, TGE_P, TDA_P, TU_P, TP_P, TK_P, 
    TG_P, TR_P, IU, IP, IK, IPD, IS, IDL, IDE, stackabilitycodedict, nbtrucks, TE, 
    TL, TW, TH, TKE, TGE, TDA, TU, TP, TK, TG, TR, TID, reverse_truckdict,  
    costinventory, costtransportation, costextratruck, timelimit
    )
    model = Model(optimizer)
    pb = TSIProblem(model, Dict{Symbol, Any}())

     
    pb[:costinventory] = costinventory
    pb[:costtransportation] = costtransportation
    pb[:costextratruck] = costextratruck
    pb[:timelimit] = timelimit
    
    pb[:item_productcodes] = item_productcodes
    pb[:truckdict] = truckdict
    pb[:supplierdict] = supplierdict
    pb[:supplierdockdict] = supplierdockdict
    pb[:plantdict] = plantdict
    pb[:plantdockdict] = plantdockdict
    pb[:nbplannedtrucks] = nbplannedtrucks
    pb[:nbitems] = nbitems
    pb[:nbsuppliers] = nbsuppliers
    pb[:nbsupplierdocks] = nbsupplierdocks
    pb[:nbplants] = nbplants
    pb[:nbplantdocks] = nbplantdocks
    # pb[:nbstacks] = nbstacks
    pb[:TE_P] = TE_P
    pb[:TL_P] = TL_P
    pb[:TW_P] = TW_P
    pb[:TH_P] = TH_P
    pb[:TKE_P] = TKE_P
    pb[:TGE_P] = TGE_P
    pb[:TDA_P] = TDA_P
    pb[:TU_P] = TU_P
    pb[:TP_P] = TP_P
    pb[:TK_P] = TK_P
    pb[:TG_P] = TG_P
    pb[:TR_P] = TR_P
    pb[:IU] = IU
    pb[:IP] = IP
    pb[:IK] = IK
    pb[:IPD] = IPD
    pb[:IS] = IS
    pb[:IDL] = IDL
    pb[:IDE] = IDE
    pb[:stackabilitycodedict] = stackabilitycodedict
    pb[:nbtrucks] = nbtrucks
    pb[:TE] = TE
    pb[:TL] = TL
    pb[:TW] = TW
    pb[:TH] = TH
    pb[:TKE] = TKE
    pb[:TGE] = TGE
    pb[:TDA] = TDA
    pb[:TU] = TU
    pb[:TP] = TP
    pb[:TK] = TK
    pb[:TG] = TG
    pb[:TR] = TR
    pb[:TID] = TID
    pb[:reverse_truckdict] = reverse_truckdict
    return pb

end

function Subproblem(t, problem, optimizer)
    
    ## Create model
    submodel = Model(optimizer)
    
    nbstacks = sum(problem[:TR][t, :])
    subproblem = Subproblem(t, problem, submodel, optimizer, nbstacks)

    # make stacks TODO

    ## Add variables
    @info "Creating variables..."
    @info "Adding zetaT..."
    @variable(submodel, zetaT >= -1)
    @info "Adding zetaE..."
    @variable(submodel, zetaE >= -1)

    # @info "Adding SS..."
    # @variable(submodel, SS[1:nbstacks] >= 0)
    @info "Adding SP..."
    @variable(submodel, SP[1:nbstacks] >= 0)
    @info "Adding SK..."
    @variable(submodel, SK[1:nbstacks] >= 0)
    @info "Adding SPD..."
    @variable(submodel, SPD[1:nbstacks] >= 0)
    @info "Adding SU..."
    @variable(submodel, SU[1:nbstacks] >= 0)
    @info "Adding SO..."
    @variable(submodel, SO[1:nbstacks] >= 0)
    @info "Adding SXe..."
    @variable(submodel, SXe[1:nbstacks] >= 0)
    @info "Adding SXo..."
    @variable(submodel, SXo[1:nbstacks] >= 0)
    @info "Adding SYe..."
    @variable(submodel, SYe[1:nbstacks] >= 0)
    @info "Adding SYo..."
    @variable(submodel, SYo[1:nbstacks] >= 0)
    @info "Adding SZe..."
    @variable(submodel, SZe[1:nbstacks] >= 0)
    @info "Adding betaM..."
    @variable(submodel, betaM[1:nbstacks] >= 0)
    @info "Adding betaP..."
    @variable(submodel, betaP[1:nbstacks] >= 0)
    @info "Adding nu..."
    @variable(submodel, nu[1:nbstacks] >= 0)
    @info "Adding tau..."
    @variable(submodel, tau[1:nbstacks] >= 0)
    @info "Adding phi..."
    @variable(submodel, phi[1:nbstacks] >= 0)
    @info "Adding SG..."
    @variable(submodel, SG[1:nbstacks, 1:nbplants] >= 0)

    @info "Adding TI..."
    @variable(submodel, TI[1:nbtrucks, 1:nbitems], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding R..."
    @variable(submodel, R[1:nbitems, 1:nbsuppliers], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding Theta..."
    @variable(submodel, Theta[1:nbitems, 1:nbsuppliers], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding S..."
    @variable(submodel, S[1:nbstacks, 1:nbitems], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding Z..."
    @variable(submodel, Z[1:nbstacks, 1:nbitems], lower_bound = 0, upper_bound = 1, Bin)

    # Omega is too big to handle, even in small instances TODO
    @debug nbstacks nbstacks
    @debug nbtrucks nbtrucks
    @debug nbitems nbitems
    @debug "nbstacks * nbtrucks * nbitems" nbstacks * nbtrucks * nbitems
    # @info "Adding Omega..."
    # @variable(submodel, Omega[1:nbstacks, 1:nbtrucks, 1:nbitems], lower_bound = 0, upper_bound = 1, container=Array, Bin)

    # @info string("Adding Omega[", t, "]...")
    # submodel[Symbol("Omega[", t, "]")] = @variable(submodel, [1:nbstacks, 1:nbitems], lower_bound = 0, upper_bound = 1, container=Array, Bin)

    # @info "Adding ST..."
    # @variable(submodel, ST[1:nbstacks, 1:nbtrucks], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding IOV..."
    @variable(submodel, IOV[1:nbitems], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding mu..."
    @variable(submodel, mu[1:nbstacks], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding eta..."
    @variable(submodel, eta[1:nbstacks], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding xi..."
    @variable(submodel, xi[1:nbstacks], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding chi..."
    @variable(submodel, chi[1:nbstacks], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding r..."
    @variable(submodel, r[1:nbstacks], lower_bound = 0, upper_bound = 1, Bin)

    @info "Adding sigma1..."
    @variable(submodel, sigma1[1:nbstacks], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding sigma2..."
    @variable(submodel, sigma2[1:nbstacks], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding sigma3..."
    @variable(submodel, sigma3[1:nbstacks], lower_bound = 0, upper_bound = 1, Bin)

    @info "Adding Psi..."
    @variable(submodel, Psi[1:nbstacks, 1:nbtrucks], lower_bound = 0, Int)
    @info "Adding Q..."
    @variable(submodel, Q[1:nbstacks, 1:nbitems], lower_bound = 0, Int)
    @info "Adding H..."
    @variable(submodel, H[1:nbstacks, 1:nbitems], lower_bound = 0, Int)
    @info "Adding V..."
    @variable(submodel, V[1:nbstacks, 1:nbitems], lower_bound = 0, Int)
    @info "Adding W..."
    @variable(submodel, W[1:nbstacks, 1:nbitems], lower_bound = 0, Int)
    @info "Adding Gl..."
    @variable(submodel, Gl[1:nbstacks, 1:nbitems], lower_bound = 0, Int)
    @info "Adding Gr..."
    @variable(submodel, Gr[1:nbstacks, 1:nbitems], lower_bound = 0, Int)

    @info "Adding lambda..."
    @variable(submodel, lambda[[1:nbstacks * (nbstacks+1)/2]], lower_bound = 0, Bin)
    
    @info "Computing parameters..."
    # MI4 = [[max(TU[:, j]...) for j in 1:first(size(TU[1, :]))] for j in 1:first(size(TI[:, 1]))]
    # MI4 = map(Int, (max(TU[:, i]...) for i = 1:first(size(TU[1, :])), j = 1:first(size(TI[:, 1]))))

    # MZ = max(IS...) + 1.0

    # MQ = max(IU...) + 1.0

    # MH = max(IP...) + 1.0

    # MV = max(IK...) + 1.0

    # MW = max(IPD...) + 1.0

    Gr = max(SO...) + 1.0 # SO? TODO

    Gl = max(IO...) + 1.0

    # Meta = nbtrucks * Mtau/10

    MTE = max(TE...)

    MTKE = max(TKE...)

    MTGE = max(TGE...)

    MTW = max(TW...) + 1.0

    # MPsi = nbitems

    Mlambda = 2.0 * TL .+ 1.0 # TODO check if Mlambda big enough

    # SZo = 0.0

    # MST = 2.0

    # MOmega = 2.0

    Mmu = 2.0

    # Mtau = 10.0

    Xi1 = falses(convert(Int, nbstacks*(nbstacks+1)/2), nbstacks)
    Xi2 = falses(convert(Int, nbstacks*(nbstacks+1)/2), nbstacks)

    epsilon = 0.001


    fillXi1!(Xi1)
    fillXi2!(Xi2)

    ## Add constraints
    @info "Adding constraints..."
    
    @info "Adding cZetaT2..."
    @constraint(submodel, cZetaT2, -zetaT >= -TIT[t] * ones(size(zetaT)[1]))

    @info "Adding cZetaE2..."
    @constraint(submodel, cZetaE2, -zetaT >= -TIE[t] * ones(size(zetaT)[1]))

    @info "Adding cTI_TR..."
    @constraint(submodel, cTI_TR, TI <= TR)

    @info "Adding cTI_1_1..."
    @constraint(submodel, cTI_1_1, transpose(TI) * ones(size(transpose(TI))[1]) == ones(size(transpose(TI))[1]))


    # @info "Adding cR_Theta_MI4..."
    # @constraint(submodel, cR_Theta_MI4, R <= Theta * MI4)

    # @info "Adding cR_TI_TU..."
    # @constraint(submodel, cR_TI_TU, -MI4*(1-Theta) <= R - (transpose(TI) * TU) <= MI4* (1-Theta))

    # @info "Adding cR_1_IU..."
    # @constraint(submodel, cR_1_IU, R * ones(size(R)[1]) == IU)

    @info "Adding cTheta_1_1..."
    @constraint(submodel, cTheta_1_1, Theta * ones(size(Theta)[1]) >= ones(size(Theta)[1]))

    # @info "Adding cIDE_TI_TDA_IDL..."
    # @constraint(submodel, cIDE_TI_TDA_IDL, IDE <= transpose(TI) * TDA <= IDL)

    # @info "Adding cZ_S_MZ..."
    # @constraint(submodel, cZ_S_MZ, Z <= S * MZ)

    # @info "Adding cPsi_Omega..."
    # @constraint(submodel, cPsi_Omega, Psi == hcat([Omega[:,i,:]*ones(size(Omega)[1]) for i in 1:size(Omega)[2]]...))

    # for j in 1:size(Omega)[2]
    #     @info "Adding cOmega_S_MOmega..."
    #     # @constraint(submodel, cOmega_S_MOmega, Omega[:, j, :] <= S * MOmega)
    #     @constraint(submodel, cOmega_S_MOmega, submodel(Symbol("Omega[", j, "]")) <= S * MOmega)
    # end
    # for i in 1:size(Omega)[1]
    #     @info "Adding cOmega_TI_S..."
    #     @constraint(submodel, cOmega_TI_S, -(1-S)*MOmega <= Omega[i,:,:] - transpose(TI) <= (1 - S)*MOmega)
    # end

    # @info "Adding cPsi_ST_MPsi..."
    # @constraint(submodel, cPsi_ST_MPsi, Psi <= ST * MPsi)

    # @info "Adding cPsi_S_ST..."
    # @constraint(submodel, cPsi_S_ST, -(1-ST)*MPsi <= Psi - hcat([S * ones(size(S)[1]) for i in 1:size(Psi)[2]]) <= (1-ST)*MPsi)

    # @info "Adding cS_IS_Z..."
    # @constraint(submodel, cS_IS_Z, S * IS == Z * ones(size(Z)[1]))

    # @info "Adding cZ_SS_S..."
    # @constraint(submodel, cZ_SS_S, -MZ * (1-S) <= Z - hcat([SS for i in 1:size(Z)[2]]) <= MZ * (1-S))

    # @info "Adding cQ_S..."
    # @constraint(submodel, cQ_S, Q <= S * MQ)

    # @info "Adding cS_IU_Q..."
    # @constraint(submodel, cS_IU_Q, S * IU == Q * ones(size(Q)[1]))

    # @info "Adding cQ_SU_S..."
    # @constraint(submodel, cQ_SU_S, -MQ * (1-S) <= Q - hcat([SU for i in 1:size(Q)[2]]) <= MQ * (1-S))


    # @info "Adding cH_S..."
    # @constraint(submodel, cH_S, H <= S * MH)

    # @info "Adding cS_IP_H..."
    # @constraint(submodel, cS_IP_H, S * IP == H * ones(size(H)[1]))

    # @info "Adding cQ_SP_S..."
    # @constraint(submodel, cQ_SP_S, -MH * (1-S) <= H - hcat([SP for i in 1:size(H)[2]]) <= MH * (1-S))


    # @info "Adding cV_S..."
    # @constraint(submodel, cV_S, V <= S * MV)

    # @info "Adding cS_IK_V..."
    # @constraint(submodel, cS_IK_V, S * IK == V * ones(size(V)[1]))

    # @info "Adding cV_SK_S..."
    # @constraint(submodel, cV_SK_S, -MV * (1-S) <= V - hcat([SK for i in 1:size(V)[2]]) <= MV * (1-S))


    # @info "Adding cW_S..."
    # @constraint(submodel, cW_S, W <= S * MW)

    # @info "Adding cS_IPD_W..."
    # @constraint(submodel, cS_IPD_W, S * IPD == W * ones(size(W)[1]))

    # @info "Adding cW_SPD_S..."
    # @constraint(submodel, cW_SPD_S, -MW * (1-S) <= W - hcat([SPD for i in 1:size(W)[2]]) <= MW * (1-S))


    @info "Adding cGl_S..."
    @constraint(submodel, cGl_S, Gl <= S * MG)
    @info "Adding cGr_S..."
    @constraint(submodel, cGr_S, Gr <= S * MG)

    @info "Adding Gl..."
    @constraint(submodel, Gl * ones(size(Gl)[1]) == Gr * ones(size(Gr)[1]))

    @info "Adding cGr_SO_S..."
    @constraint(submodel, cGr_SO_S, -MG * (1-S) <= Gr - hcat([SO for i in 1:size(Gr)[2]]) <= MG * (1-S))
    @info "Adding cGl_IOV_S..."
    @constraint(submodel, cGl_IOV_S, -MG * (1-S) <= Gl - hcat([IOV for i in 1:size(Gl)[2]]) <= MG * (1-S))

    @info "Adding cSXe_SXo_SL_SO..."
    @constraint(submodel, cSXe_SXo_SL_SO, SXe - SXo == SL + SO * MTL)
    @info "Adding cSYe_SYo_SW_SO..."
    @constraint(submodel, cSYe_SYo_SW_SO, SYe - SYo == SW + SO * MTW)

    @info "Adding cSXe_SXo_SW_SO..."
    @constraint(submodel, cSXe_SXo_SW_SO, SXe - SXo == SW + (1 - SO) * MTW)
    @info "Adding cSYe_SYo_SL_SO..."
    @constraint(submodel, cSYe_SYo_SL_SO, SYe - SYo == SL + (1 - SO) * MTL)

    @info "Adding cSZe_S_IH..."
    @constraint(submodel, cSZe_S_IH, SZe == S* IH)

    @info "Adding cSXe_ST_TL..."
    @constraint(submodel, cSXe_ST_TL, SXe <= TL[t])
    @info "Adding cSYe_ST_TW..."
    @constraint(submodel, cSYe_ST_TW, SYe <= TW[t])
    @info "Adding cSZe_ST_TH..."
    @constraint(submodel, cSZe_ST_TH, SZe <= TH[t])

    @info "Adding cSXo_SXo..."
    @constraint(submodel, cSXo_SXo, vcat(hcat([1], falses(1, size(SXo)[1]-1)), I(size(SXo)[1])) * SXo <= SXo)

    @info "Adding cXi2SXo_Xi1SXe_betaM_betaP..."
    @constraint(submodel, cXi2SXo_Xi1SXe_betaM_betaP, Xi2 * SXo - Xi1 * SXe - betaM + betaP == -epsilon)

    @info "Adding cbetaM_lambda..."
    @constraint(submodel, cbetaM_lambda, betaM <= lambda * Mlambda)

    @info "Adding betaP..."
    @constraint(submodel, betaP <= (1-lambda)*Mlambda)

    @info "Adding cmu_betaM..."
    @constraint(submodel, cmu_betaM, (1-mu) <= betaM * Mmu)

    # @info "Adding cXi2ST_Xi1ST..."
    # @constraint(submodel, cXi2ST_Xi1ST, (Xi2 * ST - Xi1 * ST) * diagm([i for i in 1:size(nu)[2]]) == nu)

    # @info "Adding ctau_phi..."
    # @constraint(submodel, ctau_phi, tau - phi <= (nu - nbtrucks) * Mtau)

    # @info "Adding ctau..."
    # @constraint(submodel, ctau, tau >= (nu-nbtrucks)*Mtau/10)

    # @info "Adding ctau_eta..."
    # @constraint(submodel, ctau_eta, tau <= eta * Meta)

    # @info "Adding cphi_eta..."
    # @constraint(submodel, cphi_eta, phi <= (1 - eta)*Meta)

    @info "Adding cXi1SYe_Xi2SYo..."
    @constraint(submodel, cXi1SYe_Xi2SYo, Xi1 * SYe <= Xi2 * SYo + xi * MTW + (1 - mu) * MTW)
    @info "Adding cXi2SYe_Xi1SYo..."
    @constraint(submodel, cXi2SYe_Xi1SYo, Xi2 * SYe <= Xi1 * SYo + (1-xi) * MTW + (1 - mu) * MTW)

    @info "Adding cXi1SU_Xi2SU..."
    @constraint(submodel, cXi1SU_Xi2SU, Xi1 * SU * TE <= Xi2 * SU * TE)

    @info "Adding cXi1SU_Xi2SU_chi..."
    @constraint(submodel, cXi1SU_Xi2SU_chi, Xi1SU - Xi2SU >= chi * epsilon - r*MTE - (1 - sigma1) * MTE)

    @info "Adding cXi2SU_Xi1SU_chi..."
    @constraint(submodel, cXi2SU_Xi1SU_chi, Xi2SU - Xi1SU >= (1-chi) * epsilon - r*MTE - (1 - sigma1) * MTE)

    @info "Adding cXi2SK_Xi1SK..."
    @constraint(submodel, cXi2SK_Xi1SK, Xi2*SK*TKE >= Xi1*SK*TKE - (1 - r) * MTKE)

    @info "Adding cXi1SK_Xi2SK_chi..."
    @constraint(submodel, cXi1SK_Xi2SK_chi, Xi1*SK*TKE - Xi2*SK*TKE >= chi*epsilon - (1 - sigma2)*MTKE)
    @info "Adding cXi2SK_Xi1SK_chi..."
    @constraint(submodel, cXi2SK_Xi1SK_chi, Xi2*SK*TKE - Xi1*SK*TKE >= (1-chi)*epsilon - (1 - sigma2)*MTKE)

    @info "Adding cXi2SG_Xi1SG..."
    @constraint(submodel, cXi2SG_Xi1SG, Xi2*SG*TGE >= Xi1*SG*TGE - (1 - sigma3) * MTGE)

    @info "Adding csigma1_sigma2_sigma3..."
    @constraint(submodel, csigma1_sigma2_sigma3, sigma1 + sigma2 + sigma3 >= 1)

    # costinventory
    # costtransportation
    # costextratruck
    # timelimit

    @objective(submodel, Min, costtransportation * zetaT + costextratruck * zetaE + costinventory * (IDL - transpose(TI[t]) * TDA) + (kappa[t] - sum([kappa[t] for t in 1:nbtrucks]) * TI)) # TODO find cost coefficients

    return subproblem
end