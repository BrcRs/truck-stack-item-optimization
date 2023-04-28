using JuMP
using MathOptInterface

using Base.Threads

include("instance_loader.jl")
include("matrix_ops.jl")


mutable struct TSIProblem
    opt::Any # Optimizer
    obj_dict::Dict{Symbol,Any}
    
end

mutable struct Subproblem
    t::Integer
    problem::TSIProblem
    submodel::Model
    optimizer::Any # Optimizer
    nbstacks::Integer
    valueTI::Matrix{Union{Missing, Bool}}
end

object_dictionary(pb::TSIProblem) = pb.obj_dict

optimizer(pb::TSIProblem) = pb.opt

model(sub::Subproblem) = sub.submodel

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

function upd_penalization!(subpb::Subproblem, TIbar, kappa)
    submodel = model(subpb)
    @objective(submodel, Min, 
    subpb[:costtransportation] * submodel[:zetaT] + 
    subpb[:costextratruck] * submodel[:zetaE] + 
    subpb[:costinventory] * (subpb[:IDL] - transpose(submodel[:TI][t]) * subpb[:TDA]) + 
    kappa * (submodel[:TI] - TIbar))
end

valueTI(subpb::Subproblem) = subpb.valueTI

function setvalueTI!(subpb::Subproblem, value::BitMatrix)
    subpb.valueTI .= value
end

function upd_valueTI!(subpb::Subproblem)
    setvalueTI!(subpb, value.(model(subpb)[:TI]))
end

function are_all_TI_equal(TIvalues, TIbar, nbtrucks, eps)
    @info "Verifying equality..."
    TIequality = true
    for t in 1:nbtrucks
        @debug "t" t
        @debug "sum(TIvalues[t, :, :] .- TIbar)" sum(abs.(TIvalues[t, :, :] .- TIbar))
        if sum(abs.(TIvalues[t, :, :] .- TIbar)) > eps
            TIequality = false
            @debug "FALSE"
            break
        end
    end
    return TIequality
end

function solve_uzawa!(problem::TSIProblem, delta::Real, eps, batchsize)
    
    nbtrucks = problem[:nbtrucks]
    nbitems = problem[:nbitems]
    # kappas = ones(nbtrucks, nbtrucks, nbitems)
    kappas = ones(nbtrucks)
    # 1. Make first solution by distributing items in planned trucks
    # Allocating TIbar
    # TIbar = vcat(problem[:TR_P], falses(problem[:nbtrucks] - problem[:nbplannedtrucks], problem[:nbitems]))
    TIvalues = falses(nbtrucks, nbtrucks, nbitems)
    TIbar = falses(nbtrucks, nbitems)
    # For each item, assign a random allowed truck
    for i in 1:nbitems
        icandidates = findall((x) -> x == 1, problem[:TR_P][:, i])
        TIbar[rand(icandidates)] = 1
    end

    # 2. Instanciate as many subproblems than there are trucks. 
    # In each subproblem, only create constraints related to the corresponding truck.
    subproblems = [Subproblem(t, problem, optimizer(problem)) for t in 1:batchsize]

    
    TIequality = are_all_TI_equal(TIvalues, TIbar, nbtrucks, eps)

    while !TIequality
        @time begin
            # @info "Convergence gap" sum(abs.(TIvalues[t, :, :] .- TIbar))
            t = 1
            while t <= nbtrucks
                for b in 0:batchsize-1
                    tb = t + b
                    @info "truck" tb
                    # If the conscensus is that there are no items in this truck, the 
                    # truck agrees and there is no need to solve it.
                    if sum(TIbar[tb, :]) == 0
                        @info "=> Truck empty"
                        # setvalueTI!(subproblems[t], TIbar)
                        TIvalues[tb, :, :] .= TIbar
                    else
                        @info "=> Optimizing subproblem..."
                        if truck(subproblems[b]) != tb
                            changetruck!(tb, subproblems[b])
                        end
                        # Solve the deterministic minimization problem for truck t with
                        # a penalization + kappa[t, k] * (TI[t, k+1] - TIbar[k])
                        # and obtain optimal first decision TI[t, k+1]
                        @info "Adding penalization..."
                        upd_penalization!(subproblems[b], TIbar, kappas[tb])
                        @info "Setting start values..."
                        set_start_value.(model(subproblems[b])[:TI], TIbar)
                        @info "Optimizing..."
                        optimize!(model(subproblems[b]))
                        # upd_valueTI!(subproblems[t])
                        TIvalues[tb, :, :] .= value.(model(subproblems[b])[:TI])
                    end
                end
                t += batchsize
            end
            # Update the mean first decisions:
            # TIbar[k+1] = sum([pi[t] * TI[t, k+1] for t in bold_T])
            # TIbar = sum([value.(model(subproblems[t])[:TI]) for t in 1:nbtrucks])
            TIbar = sum([TIvalues[t, :, :] for t in 1:nbtrucks])

            # Update the multipliers by
            for t in 1:nbtrucks
                kappas[t] = kappas[t] .+ delta * sum(TIvalues[t, :, :] .- TIbar)
            end
            TIequality = are_all_TI_equal(TIvalues, TIbar, nbtrucks, eps)
        end
    end
end

function TSIProblem(optimizer, instancepath::String)
    item_productcodes, truckdict, supplierdict, supplierdockdict, plantdict, 
    plantdockdict, nbplannedtrucks, nbitems, nbsuppliers, nbsupplierdocks, nbplants, 
    nbplantdocks, TE_P, TL_P, TW_P, TH_P, TKE_P, TGE_P, TDA_P, TU_P, TP_P, TK_P, 
    TG_P, TR_P, IU, IP, IK, IPD, IS, _IO, IL, IW, IH, IDL, IDE, stackabilitycodedict, nbtrucks, TE, 
    TL, TW, TH, TKE, TGE, TDA, TU, TP, TK, TG, TR, TID, reverse_truckdict,
    costinventory, costtransportation, costextratruck, timelimit = loadinstance(instancepath)
    

    return TSIProblem(optimizer, 
    item_productcodes, truckdict, supplierdict, supplierdockdict, plantdict, 
    plantdockdict, nbplannedtrucks, nbitems, nbsuppliers, nbsupplierdocks, nbplants, 
    nbplantdocks, TE_P, TL_P, TW_P, TH_P, TKE_P, TGE_P, TDA_P, TU_P, TP_P, TK_P, 
    TG_P, TR_P, IU, IP, IK, IPD, IS, _IO, IL, IW, IH, IDL, IDE, stackabilitycodedict, nbtrucks, TE, 
    TL, TW, TH, TKE, TGE, TDA, TU, TP, TK, TG, TR, TID, reverse_truckdict,  
    costinventory, costtransportation, costextratruck, timelimit
    )

end

getmodel(pb::TSIProblem) = pb.model

function TSIProblem(optimizer, 
    item_productcodes, truckdict, supplierdict, supplierdockdict, plantdict, 
    plantdockdict, nbplannedtrucks, nbitems, nbsuppliers, nbsupplierdocks, nbplants, 
    nbplantdocks, TE_P, TL_P, TW_P, TH_P, TKE_P, TGE_P, TDA_P, TU_P, TP_P, TK_P, 
    TG_P, TR_P, IU, IP, IK, IPD, IS, _IO, IL, IW, IH, IDL, IDE, stackabilitycodedict, nbtrucks, TE, 
    TL, TW, TH, TKE, TGE, TDA, TU, TP, TK, TG, TR, TID, reverse_truckdict,  
    costinventory, costtransportation, costextratruck, timelimit
    )
    pb = TSIProblem(optimizer, Dict{Symbol, Any}())

     
    pb[:costinventory] = costinventory
    pb[:costtransportation] = costtransportation
    pb[:costextratruck] = costextratruck
    pb[:timelimit] = timelimit
    pb[:_IO] = _IO
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
    pb[:IL] = IL
    pb[:IW] = IW
    pb[:IH] = IH
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


function changetruck!(t, subproblem::Subproblem)
    
    ## Create model
    submodel = model(subproblem)
    nbcandidateitems = sum(problem[:TR][t, :])
    # nbcandidateitems = max([sum(problem[:TR][t2, :]) for t2 in 1:nbtrucks]...)
    nbstacks = nbcandidateitems
    nbitems = problem[:nbitems]
    nbtrucks = problem[:nbtrucks]
    nbplannedtrucks = problem[:nbplannedtrucks]
    # nbplants = problem[:nbplants]
    nbplantdocks = problem[:nbplantdocks]
    nbsuppliers = problem[:nbsuppliers]
    nbsupplierdocks = problem[:nbsupplierdocks]
    nbextratrucks = nbtrucks - nbplannedtrucks
    subproblem.t = t
    # subproblem = Subproblem(t, problem, submodel, optimizer, nbstacks, Matrix{Union{Missing, Bool}}(missing, nbtrucks, nbtrucks))


    @info "Replacing S..."
    delete(submodel, S)
    unregister(submodel, S)
    @variable(submodel, S[1:nbstacks, 1:nbcandidateitems], lower_bound = 0, upper_bound = 1, Bin)

    @info "Computing parameters..."
    # MI4 = [[max(TU[:, j]...) for j in 1:first(size(TU[1, :]))] for j in 1:first(size(TI[:, 1]))]
    # MI4 = map(Int, (max(TU[:, i]...) for i = 1:first(size(TU[1, :])), j = 1:first(size(TI[:, 1]))))

    MZ = max(problem[:IS]...) + 1.0

    MQ = max(problem[:IU]...) + 1.0

    # MH = max(IP...) + 1.0

    MV = max(problem[:IK]...) + 1.0

    MW = max(problem[:IPD]...) + 1.0

    MG = max(skipmissing(problem[:_IO])...) + 1.0

    MDL =  max(problem[:IL]...) + 1.0
    MDW =  max(problem[:IW]...) + 1.0

    # Meta = nbtrucks * Mtau/10

    # MTL = Matrix{Float64}(undef, nbtrucks, 1)
    # MTW = Matrix{Float64}(undef, nbtrucks, 1)

    MTL = max(problem[:TL]...) + 1.0
    MTW = max(problem[:TW]...) + 1.0
    # if true
    # # if typeof(MTL) != Vector{T} where {T <: Real}
    #     throw(TypeError(MTL, "MTL must be of type Vector{T} where {T <: Real}", Vector{T} where {T <: Real}, typeof(MTL)))
    # end

    MTE = max(skipmissing(problem[:TE])...)

    MTKE = max(problem[:TKE]...)

    MTGE = max(problem[:TGE]...)

    MTW = max(problem[:TW]...) + 1.0

    # MPsi = nbitems

    # Mlambda = 2.0 * problem[:TL][t] + 1.0
    Mlambda = 2.0 * max(problem[:TL]...) + 1.0


    # SZo = 0.0

    # MST = 2.0

    # MOmega = 2.0

    Mmu = 2.0

    # Mtau = 10.0

    MS = 2.0

    Xi1 = falses(convert(Int, nbstacks*(nbstacks+1)/2), nbstacks)
    Xi2 = falses(convert(Int, nbstacks*(nbstacks+1)/2), nbstacks)

    epsilon = 0.001


    fillXi1!(Xi1)
    fillXi2!(Xi2)

    ## Add constraints
    @info "Replacing constraints..."
    if t <= nbplannedtrucks
        @info "Replacing cZetaT2..."
        # with_logger(ConsoleLogger(stdout, Logging.Debug, show_limited=true)) do
        #     @debug "-reshape(TI[t, :], nbitems, 1)" -reshape(TI[t, :], nbitems, 1)
        #     @debug "vones(Int8, nbitems)" vones(Int8, nbitems)
        # end
        delete(submodel, cZetaT2)
        unregister(submodel, cZetaT2)
        @constraint(submodel, cZetaT2, -zetaT .>= -TI[1:nbplannedtrucks, :] * vones(Int8, nbitems))
    else
        @info "Replacing cZetaE2..."
        delete(submodel, cZetaE2)
        unregister(submodel, cZetaE2)
        @constraint(submodel, cZetaE2, -zetaE .>= -TI[nbplannedtrucks+1:end, :] * vones(Int8, nbitems))
    end
    @info "Replacing cTI_TR..."
    delete(submodel, cTI_TR)
    unregister(submodel, cTI_TR)
    @constraint(submodel, cTI_TR, TI[t, :] <= problem[:TR][t, :])

    icandidates = findall((x) -> x == 1, problem[:TR][t, :])
    @info "Replacing cTI_1_1..."
    delete(submodel, cTI_1_1)
    unregister(submodel, cTI_1_1)
    @constraint(submodel, cTI_1_1, transpose(TI)[filter(x -> x in icandidates, 1:nbitems), :] * vones(Int8, nbtrucks) .<= vones(Int8, size(transpose(TI)[filter(x -> x in icandidates, 1:nbitems), :], 1)))
    # @debug "cTI_1_1" cTI_1_1[icandidates[1], :]
    @info "Replacing cS_TI..."
    delete(submodel, cS_TI)
    unregister(submodel, cS_TI)
    @constraint(submodel, cS_TI, S * vones(Int8, length(icandidates)) .== TI[t, filter(x -> x in icandidates, 1:nbitems)])


    @info "Replacing cS_TI_MS..."
    delete(submodel, cS_TI_MSleft)
    unregister(submodel, cS_TI_MSleft)
    @constraint(submodel, cS_TI_MSleft, -(vones(Int8, nbcandidateitems) - TI[t, filter(x -> x in icandidates, 1:nbitems)]) * MS .<= transpose(S) * vones(Int8, nbcandidateitems) - vones(Int8, nbcandidateitems))
    delete(submodel, cS_TI_MSright)
    unregister(submodel, cS_TI_MSright)
    @constraint(submodel, cS_TI_MSright, transpose(S) * vones(Int8, nbcandidateitems) - vones(Int8, nbcandidateitems) .<= (vones(Int8, nbcandidateitems) - TI[t, filter(x -> x in icandidates, 1:nbitems)]) * MS)

    @info "Replacing cS_IS_Z..."
    delete(submodel, cS_IS_Z)
    unregister(submodel, cS_IS_Z)
    @constraint(submodel, cS_IS_Z, S * problem[:IS][filter(x -> x in icandidates, 1:nbitems)] .== Z * vones(Int8, size(Z)[1]))

    @info "Replacing cS_IU_Q..."
    delete(submodel, cS_IU_Q)
    unregister(submodel, cS_IU_Q)
    @constraint(submodel, cS_IU_Q, S * problem[:IU][filter(x -> x in icandidates, 1:nbitems)] .== Q)

    @info "Replacing cS_IK_V..."
    delete(submodel, cS_IK_V)
    unregister(submodel, cS_IK_V)
    @constraint(submodel, cS_IK_V, S * problem[:IK][filter(x -> x in icandidates, 1:nbitems)] .== V)

    @info "Replacing cS_IPD_W..."
    delete(submodel, cS_IPD_W)
    unregister(submodel, cS_IPD_W)
    @constraint(submodel, cS_IPD_W, S * problem[:IPD][filter(x -> x in icandidates, 1:nbitems)] .== W)

    @info "Replacing cGl_IOV_S..."
    delete(submodel, cGl_IOV_Sleft)
    unregister(submodel, cGl_IOV_Sleft)
    @constraint(submodel, cGl_IOV_Sleft, -MG * (-S .+ 1) .<= Gl .- hcat([IOV[filter(x -> x in icandidates, 1:nbitems)] for i in 1:size(Gl)[2]]...))
    delete(submodel, cGl_IOV_Sright)
    unregister(submodel, cGl_IOV_Sright)
    @constraint(submodel, cGl_IOV_Sright, Gl .- hcat([IOV[filter(x -> x in icandidates, 1:nbitems)] for i in 1:size(Gl)[2]]...) .<= MG * (-S .+ 1))

    @info "Replacing cDL_S_IL..."
    delete(submodel, cDL_S_IL)
    unregister(submodel, cDL_S_IL)
    @constraint(submodel, cDL_S_IL, DL * vones(Int8, size(S, 1)) .== S * problem[:IL][filter(x -> x in icandidates, 1:nbitems)])
    @info "Replacing cDW_S_IW..."
    delete(submodel, cDW_S_IW)
    unregister(submodel, cDW_S_IW)
    @constraint(submodel, cDW_S_IW, DW * vones(Int8, size(S, 1)) .== S * problem[:IW][filter(x -> x in icandidates, 1:nbitems)])

    @info "Replacing cSZe_S_IH..."
    delete(submodel, cSZe_S_IH)
    unregister(submodel, cSZe_S_IH)
    @constraint(submodel, cSZe_S_IH, SZe .== S* problem[:IH][filter(x -> x in icandidates, 1:nbitems)])

    @info "Replacing cSXe_ST_TL..."
    delete(submodel, cSXe_ST_TL)
    unregister(submodel, cSXe_ST_TL)
    @constraint(submodel, cSXe_ST_TL, SXe .<= problem[:TL][t])
    @info "Replacing cSYe_ST_TW..."
    delete(submodel, cSYe_ST_TW)
    unregister(submodel, cSYe_ST_TW)
    @constraint(submodel, cSYe_ST_TW, SYe .<= problem[:TW][t])
    @info "Replacing cSZe_ST_TH..."
    delete(submodel, cSZe_ST_TH)
    unregister(submodel, cSZe_ST_TH)
    @constraint(submodel, cSZe_ST_TH, SZe .<= problem[:TH][t])

    @info "Replacing cXi1SU_Xi2SU..."
    notmissingTE = filter(x -> !ismissing(problem[:TE][t, x]), 1:nbsuppliers)
    delete(submodel, cXi1SU_Xi2SU)
    unregister(submodel, cXi1SU_Xi2SU)
    @constraint(submodel, cXi1SU_Xi2SU, Xi1 * SU[:, notmissingTE] * problem[:TE][t, notmissingTE] .<= Xi2 * SU[:, notmissingTE] * problem[:TE][t, notmissingTE])

    @info "Replacing cXi2SK_Xi1SK..."
    delete(submodel, cXi2SK_Xi1SK)
    unregister(submodel, cXi2SK_Xi1SK)
    notmissingTKE = filter(x -> !ismissing(problem[:TKE][t, x]), 1:nbsupplierdocks)
    @constraint(submodel, cXi2SK_Xi1SK, Xi2*SK[:, notmissingTKE]*problem[:TKE][t, notmissingTKE] .>= Xi1*SK[:, notmissingTKE]*problem[:TKE][t, notmissingTKE] - (-r .+ 1) * MTKE)

    @info "Replacing cXi1SK_Xi2SK_chi..."
    delete(submodel, cXi1SK_Xi2SK_chi)
    unregister(submodel, cXi1SK_Xi2SK_chi)
    @constraint(submodel, cXi1SK_Xi2SK_chi, Xi1*SK[:, notmissingTKE]*problem[:TKE][t, notmissingTKE] - Xi2*SK[:, notmissingTKE]*problem[:TKE][t, notmissingTKE] .>= chi*epsilon - (-sigma2 .+ 1)*MTKE)
    @info "Replacing cXi2SK_Xi1SK_chi..."
    delete(submodel, cXi2SK_Xi1SK_chi)
    unregister(submodel, cXi2SK_Xi1SK_chi)
    @constraint(submodel, cXi2SK_Xi1SK_chi, Xi2*SK[:, notmissingTKE]*problem[:TKE][t, notmissingTKE] - Xi1*SK[:, notmissingTKE]*problem[:TKE][t, notmissingTKE] .>= (-chi .+ 1)*epsilon - (-sigma2 .+ 1)*MTKE)

    @info "Replacing cXi2SG_Xi1SG..."
    delete(submodel, cXi2SG_Xi1SG)
    unregister(submodel, cXi2SG_Xi1SG)
    notmissingTGE = filter(x -> !ismissing(problem[:TGE][t, x]), 1:nbplantdocks)
    @constraint(submodel, cXi2SG_Xi1SG, Xi2*SG[:, notmissingTGE]*problem[:TGE][t, notmissingTGE] .>= Xi1*SG[:, notmissingTGE]*problem[:TGE][t, notmissingTGE] - (-sigma3 .+ 1) * MTGE)
end

function Subproblem(t, problem, optimizer)
    
    ## Create model
    submodel = Model(optimizer)
    nbcandidateitems = sum(problem[:TR][t, :])
    nbitems = problem[:nbitems]
    nbtrucks = problem[:nbtrucks]
    # nbcandidateitems = max([sum(problem[:TR][t2, :]) for t2 in 1:nbtrucks]...)
    nbstacks = nbcandidateitems
    nbplannedtrucks = problem[:nbplannedtrucks]
    nbplants = problem[:nbplants]
    nbplantdocks = problem[:nbplantdocks]
    nbsuppliers = problem[:nbsuppliers]
    nbsupplierdocks = problem[:nbsupplierdocks]
    nbextratrucks = nbtrucks - nbplannedtrucks
    subproblem = Subproblem(t, problem, submodel, optimizer, nbstacks, Matrix{Union{Missing, Bool}}(missing, nbtrucks, nbtrucks))

    ## Add variables
    @info "Creating variables..."
    @info "Adding zetaT..."
    @variable(submodel, zetaT[1:nbplannedtrucks] >= -1)
    @info "Adding zetaE..."
    @variable(submodel, zetaE[1:nbextratrucks] >= -1)

    @info "Adding SS..."
    @variable(submodel, SS[1:nbstacks] >= 0)
    @info "Adding SP..."
    @variable(submodel, SP[1:nbstacks], lower_bound = 0)
    @info "Adding SK..."
    @variable(submodel, SK[1:nbstacks, 1:nbsupplierdocks], lower_bound = 0, upper_bound = 1) # No need for Bin because it is constrained to be equal to their items' which are integer
    # @variable(submodel, SK[1:nbstacks], lower_bound = 0) # No need for Bin because it is constrained to be equal to their items' which are integer
    # @info "Adding SG..."
    # @variable(submodel, SG[1:nbstacks, 1:nbplantdocks], lower_bound = 0, upper_bound = 1)
    # @variable(submodel, SG[1:nbstacks], lower_bound = 0)
    @info "Adding SU..."
    @variable(submodel, SU[1:nbstacks, 1:nbsuppliers], lower_bound = 0, upper_bound = 1)
    # @variable(submodel, SU[1:nbstacks], lower_bound = 0)
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
    @variable(submodel, betaM[1:convert(Int, nbstacks*(nbstacks+1)/2)] >= 0)
    @info "Adding betaP..."
    @variable(submodel, betaP[1:convert(Int, nbstacks*(nbstacks+1)/2)] >= 0)
    @info "Adding nu..."
    @variable(submodel, nu[1:nbstacks] >= 0)
    @info "Adding tau..."
    @variable(submodel, tau[1:nbstacks] >= 0)
    @info "Adding phi..."
    @variable(submodel, phi[1:nbstacks] >= 0)
    @info "Adding SG..."
    @variable(submodel, SG[1:nbstacks, 1:nbplantdocks] >= 0, upper_bound = 1)

    @variable(submodel, SL[1:nbstacks] >= 0)
    @variable(submodel, SW[1:nbstacks] >= 0)

    @info "Adding TI..."
    @variable(submodel, TI[1:nbtrucks, 1:nbitems], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding R..."
    @variable(submodel, R[1:nbitems, 1:nbsuppliers], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding Theta..."
    @variable(submodel, Theta[1:nbitems, 1:nbsuppliers], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding S..."
    @variable(submodel, S[1:nbstacks, 1:nbcandidateitems], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding Z..."
    @variable(submodel, Z[1:nbstacks, 1:nbcandidateitems], lower_bound = 0, upper_bound = 1, Bin)

    @debug nbstacks nbstacks
    @debug nbtrucks nbtrucks
    @debug nbitems nbitems
    # @info "Adding Omega..."
    # @variable(submodel, Omega[1:nbstacks, 1:nbtrucks, 1:nbitems], lower_bound = 0, upper_bound = 1, container=Array, Bin)

    # @info string("Adding Omega[", t, "]...")
    # submodel[Symbol("Omega[", t, "]")] = @variable(submodel, [1:nbstacks, 1:nbitems], lower_bound = 0, upper_bound = 1, container=Array, Bin)

    # @info "Adding ST..."
    # @variable(submodel, ST[1:nbstacks, 1:nbtrucks], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding IOV..."
    @variable(submodel, IOV[1:nbitems], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding mu..."
    @variable(submodel, mu[1:convert(Int, nbstacks*(nbstacks+1)/2)], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding eta..."
    @variable(submodel, eta[1:convert(Int, nbstacks*(nbstacks+1)/2)], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding xi..."
    @variable(submodel, xi[1:convert(Int, nbstacks*(nbstacks+1)/2)], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding chi..."
    @variable(submodel, chi[1:convert(Int, nbstacks*(nbstacks+1)/2)], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding r..."
    @variable(submodel, r[1:convert(Int, nbstacks*(nbstacks+1)/2)], lower_bound = 0, upper_bound = 1, Bin)

    @info "Adding sigma1..."
    @variable(submodel, sigma1[1:convert(Int, nbstacks*(nbstacks+1)/2)], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding sigma2..."
    @variable(submodel, sigma2[1:convert(Int, nbstacks*(nbstacks+1)/2)], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding sigma3..."
    @variable(submodel, sigma3[1:convert(Int, nbstacks*(nbstacks+1)/2)], lower_bound = 0, upper_bound = 1, Bin)

    # @info "Adding Psi..."
    # @variable(submodel, Psi[1:nbstacks, 1:nbtrucks], lower_bound = 0, Int)
    @info "Adding Q..."
    @variable(submodel, Q[1:nbstacks, 1:nbsuppliers], lower_bound = 0, Int)
    # @info "Adding H..."
    # @variable(submodel, H[1:nbstacks, 1:nbitems], lower_bound = 0, Int)
    @info "Adding V..."
    @variable(submodel, V[1:nbstacks, 1:nbsupplierdocks], lower_bound = 0, Int)
    @info "Adding W..."
    @variable(submodel, W[1:nbstacks, 1:nbplantdocks], lower_bound = 0, Int)
    @info "Adding Gl..."
    @variable(submodel, Gl[1:nbstacks, 1:nbcandidateitems], lower_bound = 0, Int)
    @info "Adding Gr..."
    @variable(submodel, Gr[1:nbstacks, 1:nbcandidateitems], lower_bound = 0, Int)

    @variable(submodel, DL[1:nbstacks, 1:nbcandidateitems], lower_bound = 0)
    @variable(submodel, DW[1:nbstacks, 1:nbcandidateitems], lower_bound = 0)

    @info "Adding lambda..."
    @variable(submodel, lambda[1:convert(Int, nbstacks * (nbstacks+1)/2)], lower_bound = 0, Bin)
    
    @info "Computing parameters..."
    # MI4 = [[max(TU[:, j]...) for j in 1:first(size(TU[1, :]))] for j in 1:first(size(TI[:, 1]))]
    # MI4 = map(Int, (max(TU[:, i]...) for i = 1:first(size(TU[1, :])), j = 1:first(size(TI[:, 1]))))

    MZ = max(problem[:IS]...) + 1.0

    MQ = max(problem[:IU]...) + 1.0

    # MH = max(IP...) + 1.0

    MV = max(problem[:IK]...) + 1.0

    MW = max(problem[:IPD]...) + 1.0

    MG = max(skipmissing(problem[:_IO])...) + 1.0

    MDL =  max(problem[:IL]...) + 1.0
    MDW =  max(problem[:IW]...) + 1.0

    # Meta = nbtrucks * Mtau/10

    # MTL = Matrix{Float64}(undef, nbtrucks, 1)
    # MTW = Matrix{Float64}(undef, nbtrucks, 1)

    MTL = max(problem[:TL]...) + 1.0
    MTW = max(problem[:TW]...) + 1.0
    # if true
    # # if typeof(MTL) != Vector{T} where {T <: Real}
    #     throw(TypeError(MTL, "MTL must be of type Vector{T} where {T <: Real}", Vector{T} where {T <: Real}, typeof(MTL)))
    # end

    MTE = max(skipmissing(problem[:TE])...)

    MTKE = max(problem[:TKE]...)

    MTGE = max(problem[:TGE]...)

    MTW = max(problem[:TW]...) + 1.0

    # MPsi = nbitems

    # Mlambda = 2.0 * problem[:TL][t] + 1.0
    Mlambda = 2.0 * max(problem[:TL]...) + 1.0


    # SZo = 0.0

    # MST = 2.0

    # MOmega = 2.0

    Mmu = 2.0

    # Mtau = 10.0

    MS = 2.0

    Xi1 = falses(convert(Int, nbstacks*(nbstacks+1)/2), nbstacks)
    Xi2 = falses(convert(Int, nbstacks*(nbstacks+1)/2), nbstacks)

    epsilon = 0.001


    fillXi1!(Xi1)
    fillXi2!(Xi2)

    ## Add constraints
    @info "Adding constraints..."
    if t <= nbplannedtrucks
        @info "Adding cZetaT2..."
        # with_logger(ConsoleLogger(stdout, Logging.Debug, show_limited=true)) do
        #     @debug "-reshape(TI[t, :], nbitems, 1)" -reshape(TI[t, :], nbitems, 1)
        #     @debug "vones(Int8, nbitems)" vones(Int8, nbitems)
        # end
        @constraint(submodel, cZetaT2, -zetaT .>= -TI[1:nbplannedtrucks, :] * vones(Int8, nbitems))
    else
        @info "Adding cZetaE2..."
        @constraint(submodel, cZetaE2, -zetaE .>= -TI[nbplannedtrucks+1:end, :] * vones(Int8, nbitems))
    end
    @info "Adding cTI_TR..."
    @constraint(submodel, cTI_TR, TI[t, :] <= problem[:TR][t, :])

    icandidates = findall((x) -> x == 1, problem[:TR][t, :])
    @debug icandidates
    @info "Adding cTI_1_1..."
    # no more than one candidate item per truck
    # with_logger(ConsoleLogger(stdout, Logging.Debug, show_limited=true)) do
    #     @debug "transpose(TI)[filter(x -> x in icandidates, 1:nbitems), :]" transpose(TI)[filter(x -> x in icandidates, 1:nbitems), :]
    #     @debug "vones(Int8, nbtrucks)" vones(Int8, nbtrucks)
    # end
    @constraint(submodel, cTI_1_1, transpose(TI)[filter(x -> x in icandidates, 1:nbitems), :] * vones(Int8, nbtrucks) .<= vones(Int8, size(transpose(TI)[filter(x -> x in icandidates, 1:nbitems), :], 1)))
    # @debug "cTI_1_1" cTI_1_1[icandidates[1], :]
    @info "Adding cS_TI..."
    # with_logger(ConsoleLogger(stdout, Logging.Debug, show_limited=true)) do 
    #     @debug "S" S
    #     @debug "vones(Int8, length(icandidates))" vones(Int8, length(icandidates))
    #     @debug "S * vones(Int8, length(icandidates))" S * vones(Int8, length(icandidates))
    # end
    @constraint(submodel, cS_TI, S * vones(Int8, length(icandidates)) .== TI[t, filter(x -> x in icandidates, 1:nbitems)])
    # @debug "cS_TI" cS_TI
    # @info "Adding cR_Theta_MI4..."
    # @constraint(submodel, cR_Theta_MI4, R <= Theta * MI4)

    # @info "Adding cR_TI_TU..."
    # @constraint(submodel, cR_TI_TU, -MI4*(1-Theta) <= R - (transpose(TI) * TU) <= MI4* (1-Theta))

    # @info "Adding cR_1_IU..."
    # @constraint(submodel, cR_1_IU, R * vones(Int8, size(R)[1]) == IU)

    # @info "Adding cTheta_1_1..."
    # @constraint(submodel, cTheta_1_1, Theta * vones(Int8, size(Theta)[1]) >= vones(Int8, size(Theta)[1]))

    # @info "Adding cIDE_TI_TDA_IDL..."
    # @constraint(submodel, cIDE_TI_TDA_IDL, IDE <= transpose(TI) * TDA <= IDL)

    @info "Adding cZ_S_MZ..."
    @constraint(submodel, cZ_S_MZ, Z .<= S * MZ)
    # @info "Adding cPsi_Omega..."
    # @constraint(submodel, cPsi_Omega, Psi == hcat([Omega[:,i,:]*vones(Int8, size(Omega)[1]) for i in 1:size(Omega)[2]]...))

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
    # @constraint(submodel, cPsi_S_ST, -(1-ST)*MPsi <= Psi - hcat([S * vones(Int8, size(S)[1]) for i in 1:size(Psi)[2]]) <= (1-ST)*MPsi)

    @info "Adding cS_TI_MS..."
    @constraint(submodel, cS_TI_MSleft, -(vones(Int8, nbcandidateitems) - TI[t, filter(x -> x in icandidates, 1:nbitems)]) * MS .<= transpose(S) * vones(Int8, nbcandidateitems) - vones(Int8, nbcandidateitems))
    @constraint(submodel, cS_TI_MSright, transpose(S) * vones(Int8, nbcandidateitems) - vones(Int8, nbcandidateitems) .<= (vones(Int8, nbcandidateitems) - TI[t, filter(x -> x in icandidates, 1:nbitems)]) * MS)

    @info "Adding cS_IS_Z..."
    @constraint(submodel, cS_IS_Z, S * problem[:IS][filter(x -> x in icandidates, 1:nbitems)] .== Z * vones(Int8, size(Z)[1]))

    @info "Adding cZ_SS_S..."
    # with_logger(ConsoleLogger(stdout, Logging.Debug, show_limited=true)) do
    #     # @debug  "(-S .+ 1)" (-S .+ 1)
    #     # @debug "-MZ * (-S .+ 1)" -MZ * (-S .+ 1)
    #     @debug "Z" Z
    #     @debug "hcat([SS for i in 1:size(Z)[2]]...)" hcat([SS for i in 1:size(Z)[2]]...)
    # end
    @constraint(submodel, cZ_SS_Sleft, -MZ * (-S .+ 1) .<= Z .- hcat([SS for i in 1:size(Z)[2]]...))
    @constraint(submodel, cZ_SS_Sright, Z .- hcat([SS for i in 1:size(Z)[2]]...) .<= MZ * (-S .+ 1))

    @info "Adding cQ_S..."
    @constraint(submodel, cQ_S, Q .<= SU * MQ)

    @info "Adding cS_IU_Q..."
    @constraint(submodel, cS_IU_Q, S * problem[:IU][filter(x -> x in icandidates, 1:nbitems)] .== Q)

    @info "Adding cQ_SU_S..."
    # with_logger(ConsoleLogger(stdout, Logging.Debug, show_limited=true)) do
    #     @debug "SU" SU
    #     @debug "Q" Q
    # end
    @constraint(submodel, cQ_SU_Sleft, -MQ * (1 .- SU) .<= Q .- hcat([S * vones(Int8, size(S, 2)) for i in 1:size(Q, 2)]...))
    @constraint(submodel, cQ_SU_Sright, Q .- hcat([S * vones(Int8, size(S, 2)) for i in 1:size(Q, 2)]...) .<= MQ * (1 .- SU))


    # @info "Adding cH_S..."
    # @constraint(submodel, cH_S, H <= S * MH)

    # @info "Adding cS_IP_H..."
    # @constraint(submodel, cS_IP_H, S * IP == H * vones(Int8, size(H)[1]))

    # @info "Adding cH_SP_S..."
    # @constraint(submodel, cH_SP_S, -MH * (1-S) <= H - hcat([SP for i in 1:size(H)[2]]) <= MH * (1-S))


    @info "Adding cV_S..."
    @constraint(submodel, cV_S, V .<= SK * MV)

    @info "Adding cS_IK_V..."
    @constraint(submodel, cS_IK_V, S * problem[:IK][filter(x -> x in icandidates, 1:nbitems)] .== V)

    @info "Adding cV_SK_S..."
    @constraint(submodel, cV_SK_Sleft, -MV * (-SK .+ 1) .<= V .- hcat([S * vones(Int8, size(S, 2)) for i in 1:size(V, 2)]...))
    @constraint(submodel, cV_SK_Sright, V .- hcat([S * vones(Int8, size(S, 2)) for i in 1:size(V, 2)]...) .<= MV * (-SK .+ 1))


    @info "Adding cW_S..."
    @constraint(submodel, cW_S, W .<= SG * MW)

    @info "Adding cS_IPD_W..."

    @constraint(submodel, cS_IPD_W, S * problem[:IPD][filter(x -> x in icandidates, 1:nbitems)] .== W)

    @info "Adding cW_SG_S..."
    @constraint(submodel, cW_SPD_Sleft, -MW * (-SG .+ 1) .<= W .- hcat([S * vones(Int8, size(S, 2)) for i in 1:size(W, 2)]...))
    @constraint(submodel, cW_SPD_Sright, W .- hcat([S * vones(Int8, size(S, 2)) for i in 1:size(W, 2)]...) .<= MW * (-SG .+ 1))


    @info "Adding cGl_S..."
    @constraint(submodel, cGl_S, Gl .<= S * MG)
    @info "Adding cGr_S..."
    @constraint(submodel, cGr_S, Gr .<= S * MG)

    @info "Adding Gl..."
    @constraint(submodel, Gl * vones(Int8, size(Gl)[1]) .== Gr * vones(Int8, size(Gr)[1]))

    @info "Adding cGr_SO_S..."
    @constraint(submodel, cGr_SO_Sleft, -MG * (-S .+ 1) .<= Gr .- hcat([SO for i in 1:size(Gr)[2]]...))
    @constraint(submodel, cGr_SO_Sright, Gr .- hcat([SO for i in 1:size(Gr)[2]]...) .<= MG * (-S .+ 1))
    @info "Adding cGl_IOV_S..."
    @constraint(submodel, cGl_IOV_Sleft, -MG * (-S .+ 1) .<= Gl .- hcat([IOV[filter(x -> x in icandidates, 1:nbitems)] for i in 1:size(Gl)[2]]...))
    @constraint(submodel, cGl_IOV_Sright, Gl .- hcat([IOV[filter(x -> x in icandidates, 1:nbitems)] for i in 1:size(Gl)[2]]...) .<= MG * (-S .+ 1))

    @info "Adding cDL_S..."
    @constraint(submodel, cDL_S, DL .<= S * MDL)
    @info "Adding cDL_SL..."
    @constraint(submodel, cDL_SLleft, -MDL * (-S .+ 1) .<= DL - hcat([SL for i in 1:size(DL)[2]]...))
    @constraint(submodel, cDL_SLright, DL - hcat([SL for i in 1:size(DL)[2]]...) .<= MDL * (-S .+ 1))
    @info "Adding cDL_S_IL..."
    @constraint(submodel, cDL_S_IL, DL * vones(Int8, size(S, 1)) .== S * problem[:IL][filter(x -> x in icandidates, 1:nbitems)])
    @info "Adding cDW_S..."
    @constraint(submodel, cDW_S, DW .<= S * MDW)
    @info "Adding cDW_SW..."
    @constraint(submodel, cDW_SWleft, -MDW * (-S .+ 1) .<= DW - hcat([SW for i in 1:size(DW)[2]]...))
    @constraint(submodel, cDW_SWright, DW - hcat([SW for i in 1:size(DW)[2]]...) .<= MDW * (-S .+ 1))
    @info "Adding cDW_S_IW..."
    @constraint(submodel, cDW_S_IW, DW * vones(Int8, size(S, 1)) .== S * problem[:IW][filter(x -> x in icandidates, 1:nbitems)])

    @info "Adding cSXe_SXo_SL_SO..."
    @constraint(submodel, cSXe_SXo_SL_SO, SXe - SXo .== SL + SO * MTL)
    @info "Adding cSYe_SYo_SW_SO..."
    @constraint(submodel, cSYe_SYo_SW_SO, SYe - SYo .== SW + SO * MTW)

    @info "Adding cSXe_SXo_SW_SO..."
    @constraint(submodel, cSXe_SXo_SW_SO, SXe - SXo .== SW + (-SO .+ 1) * MTW)
    @info "Adding cSYe_SYo_SL_SO..."
    @constraint(submodel, cSYe_SYo_SL_SO, SYe - SYo .== SL + (-SO .+ 1) * MTL)

    @info "Adding cSZe_S_IH..."
    @constraint(submodel, cSZe_S_IH, SZe .== S* problem[:IH][filter(x -> x in icandidates, 1:nbitems)])

    @info "Adding cSXe_ST_TL..."
    @constraint(submodel, cSXe_ST_TL, SXe .<= problem[:TL][t])
    @info "Adding cSYe_ST_TW..."
    @constraint(submodel, cSYe_ST_TW, SYe .<= problem[:TW][t])
    @info "Adding cSZe_ST_TH..."
    @constraint(submodel, cSZe_ST_TH, SZe .<= problem[:TH][t])

    @info "Adding cSXo_SXo..."
    # with_logger(ConsoleLogger(stdout, Logging.Debug, show_limited=true)) do
    #     @debug  "SXo[1:end-1]" SXo[1:end-1]
    #     @debug "SXo[2:end]" SXo[2:end]
    # end
    # @constraint(submodel, cSXo_SXo, (vcat(hcat([1], falses(1, size(SXo)[1]-1)), I(size(SXo)[1])) * SXo)[2:end] .<= SXo)
    @constraint(submodel, cSXo_SXo, SXo[1:end-1] .<= SXo[2:end])

    @info "Adding cXi2SXo_Xi1SXe_betaM_betaP..."
    # with_logger(ConsoleLogger(stdout, Logging.Debug, show_limited=true)) do
    #     @debug "Xi2" Xi2
    #     @debug "SXo" SXo
    #     @debug "Xi1" Xi1
    #     @debug "SXe" SXe
    #     @debug "betaM" betaM
    #     @debug "betaP" betaP
    #     @debug "-epsilon * vones(Float64, size(Xi1, 1))" -epsilon * vones(Float64, size(Xi1, 1))
    # end
    @constraint(submodel, cXi2SXo_Xi1SXe_betaM_betaP, (Xi2 * SXo) - (Xi1 * SXe) - betaM + betaP .== -epsilon * vones(Float64, size(Xi1, 1)))

    @info "Adding cbetaM_lambda..."
    # with_logger(ConsoleLogger(stdout, Logging.Debug, show_limited=true)) do
    #     @debug "betaM" betaM
    #     @debug "lambda" lambda
    #     @debug "Mlambda" Mlambda 
    # end
    @constraint(submodel, cbetaM_lambda, betaM .<= lambda .* Mlambda)

    @info "Adding betaP..."
    @constraint(submodel, betaP .<= (-lambda .+ 1)*Mlambda)

    @info "Adding cmu_betaM..."
    @constraint(submodel, cmu_betaM, (-mu .+ 1) .<= betaM * Mmu)

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
    @constraint(submodel, cXi1SYe_Xi2SYo, Xi1 * SYe .<= Xi2 * SYo + xi * MTW + (-mu .+ 1) * MTW)
    @info "Adding cXi2SYe_Xi1SYo..."
    @constraint(submodel, cXi2SYe_Xi1SYo, Xi2 * SYe .<= Xi1 * SYo + (-xi .+ 1) * MTW + (-mu .+ 1) * MTW)

    @info "Adding cXi1SU_Xi2SU..."
    # with_logger(ConsoleLogger(stdout, Logging.Debug, show_limited=true)) do
    #     @debug "Xi1" Xi1
    #     @debug "SU" SU
    #     @debug "problem[:TE]" problem[:TE]
    #     @debug "problem[:TE][t, :]" problem[:TE][t, :]
    # end
    notmissingTE = filter(x -> !ismissing(problem[:TE][t, x]), 1:nbsuppliers)
    @constraint(submodel, cXi1SU_Xi2SU, Xi1 * SU[:, notmissingTE] * problem[:TE][t, notmissingTE] .<= Xi2 * SU[:, notmissingTE] * problem[:TE][t, notmissingTE])
    @info "Adding cXi1SU_Xi2SU_chi..."
    @constraint(submodel, cXi1SU_Xi2SU_chi, Xi1*SU - Xi2*SU .>= chi * epsilon - r*MTE - (-sigma1 .+ 1) * MTE)

    @info "Adding cXi2SU_Xi1SU_chi..."
    @constraint(submodel, cXi2SU_Xi1SU_chi, Xi2*SU - Xi1*SU .>= (-chi .+ 1) * epsilon - r*MTE - (-sigma1 .+ 1) * MTE)

    @info "Adding cXi2SK_Xi1SK..."
    notmissingTKE = filter(x -> !ismissing(problem[:TKE][t, x]), 1:nbsupplierdocks)
    @constraint(submodel, cXi2SK_Xi1SK, Xi2*SK[:, notmissingTKE]*problem[:TKE][t, notmissingTKE] .>= Xi1*SK[:, notmissingTKE]*problem[:TKE][t, notmissingTKE] - (-r .+ 1) * MTKE)

    @info "Adding cXi1SK_Xi2SK_chi..."
    @constraint(submodel, cXi1SK_Xi2SK_chi, Xi1*SK[:, notmissingTKE]*problem[:TKE][t, notmissingTKE] - Xi2*SK[:, notmissingTKE]*problem[:TKE][t, notmissingTKE] .>= chi*epsilon - (-sigma2 .+ 1)*MTKE)
    @info "Adding cXi2SK_Xi1SK_chi..."
    @constraint(submodel, cXi2SK_Xi1SK_chi, Xi2*SK[:, notmissingTKE]*problem[:TKE][t, notmissingTKE] - Xi1*SK[:, notmissingTKE]*problem[:TKE][t, notmissingTKE] .>= (-chi .+ 1)*epsilon - (-sigma2 .+ 1)*MTKE)

    @info "Adding cXi2SG_Xi1SG..."
    notmissingTGE = filter(x -> !ismissing(problem[:TGE][t, x]), 1:nbplantdocks)
    # with_logger(ConsoleLogger(stdout, Logging.Debug, show_limited=true)) do
    #     @debug "SG[:, notmissingTGE]" SG[:, notmissingTGE]
    #     @debug "problem[:TGE][t, notmissingTGE]" problem[:TGE][t, notmissingTGE]
    # end
    @constraint(submodel, cXi2SG_Xi1SG, Xi2*SG[:, notmissingTGE]*problem[:TGE][t, notmissingTGE] .>= Xi1*SG[:, notmissingTGE]*problem[:TGE][t, notmissingTGE] - (-sigma3 .+ 1) * MTGE)

    @info "Adding csigma1_sigma2_sigma3..."
    @constraint(submodel, csigma1_sigma2_sigma3, sigma1 + sigma2 + sigma3 .>= 1)

    @info "Adding IOV constraints..."
    for i in 1:nbitems
        if !ismissing(problem[:_IO][i])
            @constraint(submodel, IOV[i] == problem[:_IO][i])
        end
    end
    @info "Adding objective function..."
    # with_logger(ConsoleLogger(stdout, Logging.Debug, show_limited=true)) do
    #     @debug "problem[:costtransportation]" problem[:costtransportation]
    #     @debug "zetaT" zetaT

    #     @debug "transpose(TI)" transpose(TI)
    #     @debug "problem[:TDA]" problem[:TDA]
    #     @debug "problem[:costinventory]" problem[:costinventory]
    #     # @debug "problem[:IDL] - transpose(TI) * problem[:TDA]" problem[:IDL] - transpose(TI) * problem[:TDA]
    # end
    # @objective(submodel, Min, 
    # problem[:costtransportation] * zetaT + 
    # problem[:costextratruck] * zetaE + 
    # problem[:costinventory] * (problem[:IDL] - transpose(TI) * problem[:TDA]) + 
    # (kappa[t] - sum([problem[:kappa][t] for t in 1:nbtrucks]) * TI))



    @objective(submodel, Min, 
    sum(problem[:costtransportation] * zetaT) + 
    sum(problem[:costextratruck] * zetaE) + 
    sum(problem[:costinventory] * (problem[:IDL] - transpose(TI) * problem[:TDA])))

    @info "Done!"
    # @debug "kappa[t]" kappa[t]
    # @debug "sum([kappa[t2] for t2 in 1:nbtrucks])" sum([kappa[t2] for t2 in 1:nbtrucks])
    # @objective(submodel, Min, kappa * (submodel[:TI] - TIbar))
    # @info "Adding objective function 2..."
    # @objective(submodel, Min, problem[:costinventory] * (problem[:IDL] - transpose(TI) * problem[:TDA]))
    # @info "Adding objective function 3..."
    # @objective(submodel, Min, problem[:costextratruck] * zetaE + problem[:costinventory] * (problem[:IDL] - transpose(TI) * problem[:TDA]))

    return subproblem
end