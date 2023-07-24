using JuMP

include("instance_loader.jl")
include("matrix_ops.jl")

mutable struct TSIModel
    model::Model
end

function unregister(model::TSIModel, key::Symbol)
    return JuMP.unregister(model.model, key)
end

function Base.getindex(m::TSIModel, name::Symbol)
    return getindex(m.model, name)
end

function Base.setindex!(model::TSIModel, value, name::Symbol)
    return Base.setindex!(model.model, value, name)
end

function Base.haskey(model::TSIModel, name::Symbol)
    return Base.haskey(model.model, name)
end

# function Base.getindex(x::Array{<:AbstractJuMPScalar}; kwargs...)
#     return Base.getindex(x, kwargs...)
# end

function TSIModel(optimizer, instancepath)
    item_productcodes, truckdict, supplierdict, supplierdockdict, plantdict, 
    plantdockdict, nbplannedtrucks, nbitems, nbsuppliers, nbsupplierdocks, nbplants, 
    nbplantdocks, TE_P, TL_P, TW_P, TH_P, TKE_P, TGE_P, TDA_P, TU_P, TP_P, TK_P, 
    TG_P, TR_P, IU, IP, IK, IPD, IS, IDL, IDE, stackabilitycodedict, nbtrucks, TE, 
    TL, TW, TH, TKE, TGE, TDA, TU, TP, TK, TG, TR, TID, reverse_truckdict = loadinstance(instancepath)
    
    nbstacks = nbitems

    return TSIModel(optimizer, 
    item_productcodes, truckdict, supplierdict, supplierdockdict, plantdict, 
    plantdockdict, nbplannedtrucks, nbitems, nbsuppliers, nbsupplierdocks, nbplants, 
    nbplantdocks, nbstacks, TE_P, TL_P, TW_P, TH_P, TKE_P, TGE_P, TDA_P, TU_P, TP_P, TK_P, 
    TG_P, TR_P, IU, IP, IK, IPD, IS, IDL, IDE, stackabilitycodedict, nbtrucks, TE, 
    TL, TW, TH, TKE, TGE, TDA, TU, TP, TK, TG, TR, TID, reverse_truckdict
    )

end

function TSIModel(optimizer, 
    item_productcodes, truckdict, supplierdict, supplierdockdict, plantdict, 
    plantdockdict, nbplannedtrucks, nbitems, nbsuppliers, nbsupplierdocks, nbplants, 
    nbplantdocks, nbstacks, TE_P, TL_P, TW_P, TH_P, TKE_P, TGE_P, TDA_P, TU_P, TP_P, TK_P, 
    TG_P, TR_P, IU, IP, IK, IPD, IS, IDL, IDE, stackabilitycodedict, nbtrucks, TE, 
    TL, TW, TH, TKE, TGE, TDA, TU, TP, TK, TG, TR, TID, reverse_truckdict
    )

    ## Create model
    model = Model(optimizer)
    
    model[:item_productcodes] = item_productcodes
    model[:truckdict] = truckdict
    model[:supplierdict] = supplierdict
    model[:supplierdockdict] = supplierdockdict
    model[:plantdict] = plantdict
    model[:plantdockdict] = plantdockdict
    model[:nbplannedtrucks] = nbplannedtrucks
    model[:nbitems] = nbitems
    model[:nbsuppliers] = nbsuppliers
    model[:nbsupplierdocks] = nbsupplierdocks
    model[:nbplants] = nbplants
    model[:nbplantdocks] = nbplantdocks
    model[:nbstacks] = nbstacks
    model[:TE_P] = TE_P
    model[:TL_P] = TL_P
    model[:TW_P] = TW_P
    model[:TH_P] = TH_P
    model[:TKE_P] = TKE_P
    model[:TGE_P] = TGE_P
    model[:TDA_P] = TDA_P
    model[:TU_P] = TU_P
    model[:TP_P] = TP_P
    model[:TK_P] = TK_P
    model[:TG_P] = TG_P
    model[:TR_P] = TR_P
    model[:IU] = IU
    model[:IP] = IP
    model[:IK] = IK
    model[:IPD] = IPD
    model[:IS] = IS
    model[:IDL] = IDL
    model[:IDE] = IDE
    model[:stackabilitycodedict] = stackabilitycodedict
    model[:nbtrucks] = nbtrucks
    model[:TE] = TE
    model[:TL] = TL
    model[:TW] = TW
    model[:TH] = TH
    model[:TKE] = TKE
    model[:TGE] = TGE
    model[:TDA] = TDA
    model[:TU] = TU
    model[:TP] = TP
    model[:TK] = TK
    model[:TG] = TG
    model[:TR] = TR
    model[:TID] = TID
    model[:reverse_truckdict] = reverse_truckdict
    


    ## Add variables
    @info "Creating variables..."
    @info "Adding zetaT..."
    @variable(model, zetaT[1:nbplannedtrucks] >= 0)
    @info "Adding zetaE..."
    @variable(model, zetaE[1:nbtrucks - nbplannedtrucks] >= 0)

    @info "Adding SS..."
    @variable(model, SS[1:nbstacks] >= 0)
    @info "Adding SP..."
    @variable(model, SP[1:nbstacks] >= 0)
    @info "Adding SK..."
    @variable(model, SK[1:nbstacks] >= 0)
    @info "Adding SPD..."
    @variable(model, SPD[1:nbstacks] >= 0)
    @info "Adding SU..."
    @variable(model, SU[1:nbstacks] >= 0)
    @info "Adding SO..."
    @variable(model, SO[1:nbstacks] >= 0)
    @info "Adding SXe..."
    @variable(model, SXe[1:nbstacks] >= 0)
    @info "Adding SXo..."
    @variable(model, SXo[1:nbstacks] >= 0)
    @info "Adding SYe..."
    @variable(model, SYe[1:nbstacks] >= 0)
    @info "Adding SYo..."
    @variable(model, SYo[1:nbstacks] >= 0)
    @info "Adding SZe..."
    @variable(model, SZe[1:nbstacks] >= 0)
    @info "Adding betaM..."
    @variable(model, betaM[1:nbstacks] >= 0)
    @info "Adding betaP..."
    @variable(model, betaP[1:nbstacks] >= 0)
    @info "Adding nu..."
    @variable(model, nu[1:nbstacks] >= 0)
    @info "Adding tau..."
    @variable(model, tau[1:nbstacks] >= 0)
    @info "Adding phi..."
    @variable(model, phi[1:nbstacks] >= 0)
    @info "Adding SG..."
    @variable(model, SG[1:nbstacks, 1:nbplants] >= 0)

    @info "Adding TI..."
    @variable(model, TI[1:nbtrucks, 1:nbitems], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding R..."
    @variable(model, R[1:nbitems, 1:nbsuppliers], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding Theta..."
    @variable(model, Theta[1:nbitems, 1:nbsuppliers], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding S..."
    @variable(model, S[1:nbstacks, 1:nbitems], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding Z..."
    @variable(model, Z[1:nbstacks, 1:nbitems], lower_bound = 0, upper_bound = 1, Bin)

    # Omega is too big to handle, even in small instances TODO
    @debug nbstacks nbstacks
    @debug nbtrucks nbtrucks
    @debug nbitems nbitems
    @debug "nbstacks * nbtrucks * nbitems" nbstacks * nbtrucks * nbitems
    @info "Adding Omega..."
    # @variable(model, Omega[1:nbstacks, 1:nbtrucks, 1:nbitems], lower_bound = 0, upper_bound = 1, container=Array, Bin)

    # TODO this is taking way too long
    for t in 1:nbtrucks
        @info string("Adding Omega[", t, "]...")
        model[Symbol("Omega[", t, "]")] = @variable(model, [1:nbstacks, 1:nbitems], lower_bound = 0, upper_bound = 1, container=Array, Bin)
    end

    @info "Adding ST..."
    @variable(model, ST[1:nbstacks, 1:nbtrucks], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding IOV..."
    @variable(model, IOV[1:nbitems], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding mu..."
    @variable(model, mu[1:nbstacks], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding eta..."
    @variable(model, eta[1:nbstacks], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding xi..."
    @variable(model, xi[1:nbstacks], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding chi..."
    @variable(model, chi[1:nbstacks], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding r..."
    @variable(model, r[1:nbstacks], lower_bound = 0, upper_bound = 1, Bin)

    @info "Adding sigma1..."
    @variable(model, sigma1[1:nbstacks], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding sigma2..."
    @variable(model, sigma2[1:nbstacks], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding sigma3..."
    @variable(model, sigma3[1:nbstacks], lower_bound = 0, upper_bound = 1, Bin)

    @info "Adding Psi..."
    @variable(model, Psi[1:nbstacks, 1:nbtrucks], lower_bound = 0, Int)
    @info "Adding Q..."
    @variable(model, Q[1:nbstacks, 1:nbitems], lower_bound = 0, Int)
    @info "Adding H..."
    @variable(model, H[1:nbstacks, 1:nbitems], lower_bound = 0, Int)
    @info "Adding V..."
    @variable(model, V[1:nbstacks, 1:nbitems], lower_bound = 0, Int)
    @info "Adding W..."
    @variable(model, W[1:nbstacks, 1:nbitems], lower_bound = 0, Int)
    @info "Adding Gl..."
    @variable(model, Gl[1:nbstacks, 1:nbitems], lower_bound = 0, Int)
    @info "Adding Gr..."
    @variable(model, Gr[1:nbstacks, 1:nbitems], lower_bound = 0, Int)

    @info "Adding lambda..."
    @variable(model, lambda[[1:nbstacks * (nbstacks+1)/2]], lower_bound = 0, Bin)
    
    @info "Computing parameters..."
    # MI4 = [[max(TU[:, j]...) for j in 1:first(size(TU[1, :]))] for j in 1:first(size(TI[:, 1]))]
    MI4 = map(Int, (max(TU[:, i]...) for i = 1:first(size(TU[1, :])), j = 1:first(size(TI[:, 1]))))

    MZ = max(IS...) + 1.0

    MQ = max(IU...) + 1.0

    MH = max(IP...) + 1.0

    MV = max(IK...) + 1.0

    MW = max(IPD...) + 1.0

    Gr = max(SO...) + 1.0 # SO? TODO

    Gl = max(IO...) + 1.0

    Meta = nbtrucks * Mtau/10

    MTE = max(TE...)

    MTKE = max(TKE...)

    MTGE = max(TGE...)

    MTW = max(TW...) + 1.0

    MPsi = nbitems

    Mlambda = 2.0 * TL .+ 1.0 # TODO check if Mlambda big enough

    SZo = 0.0

    MST = 2.0

    MOmega = 2.0

    Mmu = 2.0

    Mtau = 10.0

    Xi1 = falses(convert(Int, nbstacks*(nbstacks+1)/2), nbstacks)
    Xi2 = falses(convert(Int, nbstacks*(nbstacks+1)/2), nbstacks)

    epsilon = 0.001


    fillXi1!(Xi1)
    fillXi2!(Xi2)

    ## Add constraints
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
    @constraint(model, cPsi_ST_MPsi, Psi <= ST * MPsi)

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
    @constraint(model, ctau_phi_nu, tau - phi <= (nu - nbtrucks) * Mtau)

    @info "Adding ctau_nu..."
    @constraint(model, ctau_nu, tau >= (nu-nbtrucks)*Mtau/10)

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
    return model
end