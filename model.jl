using JuMP
using MathOptInterface

using Base.Threads
# using Clp # Only for debug
# using CDD # Only debug
include("instance_loader.jl")
include("matrix_ops.jl")
include("linear_infeasibilities.jl")

# TODO Refactor

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

truck(sub::Subproblem) = sub.t

function JuMP.unregister(pb::TSIProblem, key::Symbol)
    return JuMP.JuMP.unregister(pb.model, key)
end

function Base.getindex(p::TSIProblem, name::Symbol)
    obj_dict = object_dictionary(p)
    if !haskey(obj_dict, name)
        throw(KeyError(name))
    end
    return obj_dict[name]
end

function Base.getindex(sub::Subproblem, name::Symbol)
    obj_dict = object_dictionary(sub.problem)
    if !haskey(obj_dict, name)
        throw(KeyError(name))
    end
    return obj_dict[name]
end


function Base.setindex!(problem::TSIProblem, value, name::Symbol)
    return object_dictionary(problem)[name] = value
end

function JuMP.unregister(problem::TSIProblem, key::Symbol)
    delete!(object_dictionary(problem), key)
    return
end
function Base.haskey(problem::TSIProblem, name::Symbol)
    return haskey(object_dictionary(problem), name)
end

function Base.haskey(subproblem::Subproblem, name::Symbol)
    return haskey(object_dictionary(subproblem.problem), name)
end

function upd_penalization!(subpb::Subproblem, TIbar, kappa)
    submodel = model(subpb)
    t = truck(subpb)

    nbplannedtrucks = subpb[:nbplannedtrucks]
    nbitems = subpb[:nbitems]
    # @objective(submodel, Min, 
    # sum(problem[:costtransportation] * zetaT) + 
    # sum(problem[:costextratruck] * zetaE) + 
    # sum(problem[:costinventory] * (problem[:IDL] - transpose(TI) * problem[:TDA])))

    # @debug "subpb[:costtransportation]" subpb[:costtransportation]
    # @debug "submodel[:zetaT]" submodel[:zetaT]
    # @debug "subpb[:costextratruck]" subpb[:costextratruck]
    # @debug "submodel[:zetaE]" submodel[:zetaE]
    # @debug "subpb[:costinventory]" subpb[:costinventory]
    # @debug "subpb[:IDL]" subpb[:IDL]
    # @debug "submodel[:TI][t, :]" submodel[:TI][t, :]
    # @debug "subpb[:TDA][t]" subpb[:TDA][t]
    # @debug "kappa" kappa
    # @debug "submodel[:TI]" submodel[:TI]
    # @debug "TIbar" TIbar
    # @debug "submodel[:TI][t, :] * subpb[:TDA][t]" submodel[:TI][t, :] * subpb[:TDA][t]
    # @debug "subpb[:IDL] - submodel[:TI][t, :] * subpb[:TDA][t]" subpb[:IDL] - submodel[:TI][t, :] * subpb[:TDA][t]
    # @debug "sum(subpb[:IDL] - submodel[:TI][t, :] * subpb[:TDA][t])" sum(subpb[:IDL] - submodel[:TI][t, :] * subpb[:TDA][t])

    # error("Stop")

    # Worst case solution value:
    MGI = subpb[:costtransportation] * nbplannedtrucks + (nbitems - nbplannedtrucks) * subpb[:costextratruck] + subpb[:costinventory] * nbitems * (min(subpb[:IDL]...) - max(subpb[:TDE]...))

    @expression(submodel, obj1, sum(subpb[:costtransportation] * submodel[:zetaT]))
    @expression(submodel, obj2, sum(subpb[:costextratruck] * submodel[:zetaE]))
    @expression(submodel, obj3, subpb[:costinventory] * sum(subpb[:IDL] - submodel[:TI][t, :] * subpb[:TDA][t]))
    @expression(submodel, obj4, kappa * sum(submodel[:TI] - TIbar))
    @expression(submodel, obj5, MGI * sum(submodel[:GI]))

    # @objective(submodel, Min, 
    # sum(subpb[:costtransportation] * submodel[:zetaT]) + 
    # sum(subpb[:costextratruck] * submodel[:zetaE]) + 
    # subpb[:costinventory] * sum(subpb[:IDL] - submodel[:TI][t, :] * subpb[:TDA][t]) + 
    # kappa * sum(submodel[:TI] - TIbar))

    @objective(submodel, Min, obj1 + obj2 + obj3 + obj4 + obj5)

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
        # @debug "t" t
        # @debug "sum(TIvalues[t, :, :] .- TIbar)" sum(abs.(TIvalues[t, :, :] .- TIbar))
        if sum(TIvalues[t, :, :] == TIbar ? 0.0 : abs.(TIvalues[t, :, :] .- TIbar)) > eps
            TIequality = false
            # @debug "FALSE"
            break
        end
    end
    return TIequality
end

function solve_uzawa!(problem::TSIProblem, delta::Real, eps, batchsize, chosentrucks)
    nbtrucks = problem[:nbtrucks]
    nbplannedtrucks = problem[:nbplannedtrucks]
    nbitems = problem[:nbitems]
    nbchosentrucks = length(chosentrucks)
    # kappas = ones(nbtrucks, nbtrucks, nbitems)
    kappas = ones(nbchosentrucks)
    # 1. Make first solution by distributing items in planned trucks
    # Allocating TIbar
    # TIbar = vcat(problem[:TR_P], falses(problem[:nbchosentrucks] - problem[:nbplannedtrucks], problem[:nbitems]))
    TIvalues = falses(nbchosentrucks, nbchosentrucks, nbitems) # TODO might need to be shared array
    TIbar = falses(nbchosentrucks, nbitems)
    # For each item, assign a random allowed truck
    for i in 1:nbitems
        tcandidates = findall((x) -> x == 1, problem[:TR_P][:, i])
        TIbar[rand(tcandidates), i] = 1
    end

    # 2. Instanciate as many subproblems than there are trucks. 
    # In each subproblem, only create constraints related to the corresponding truck.
    subproblems = [Subproblem(t, problem, optimizer(problem), chosentrucks) for t in 1:batchsize]

    optsolpertruck = [Dict{Symbol, Any}() for t in 1:nbtrucks]
    
    TIequality = are_all_TI_equal(TIvalues, TIbar, nbchosentrucks, eps)

    # alltrucksdone = falses(nbchosentrucks)
    firstpass = true # debug purpose
    while !TIequality || firstpass
        firstpass = false # debug purpose
        @time begin
            # @info "Convergence gap" sum(abs.(TIvalues[t, :, :] .- TIbar))
            pt = 1

            for pt in 1:nbplannedtrucks
                # @sync @distributed for b in 0:batchsize-1 # Threads
                for b in 0:batchsize-1 # Threads
                    @debug "Thread $b" 
                    ### Assign a planned truck (and corresponding extra trucks) to a thread/subproblem instance
                    ptb = pt + b
                    if ptb > nbplannedtrucks
                        continue
                    end
                    for i in filter(x -> x in chosentrucks, problem[:truckindices][ptb])
                        @debug "truck" i
                        # If the conscensus is that there are no items in this truck, the 
                        # truck agrees and there is no need to solve it.
                        if sum(TIbar[i, :]) == 0
                            @debug "=> Truck empty"
                            # setvalueTI!(subproblems[t], TIbar)
                            TIvalues[i, :, :] .= TIbar
                            # optsolpertruck[i][:S] = nothing
                            # optsolpertruck[i][:SG] = nothing
                            # optsolpertruck[i][:SL] = nothing
                            # optsolpertruck[i][:SW] = nothing
                            # optsolpertruck[i][:SH] = nothing

                        else
                            @debug "=> Optimizing subproblem..."
                            if truck(subproblems[b+1]) != i
                                changetruck!(i, subproblems[b+1], chosentrucks, changeS=i == ptb)
                            end
                            # Solve the deterministic minimization problem for truck t with
                            # a penalization + kappa[t, k] * (TI[t, k+1] - TIbar[k])
                            # and obtain optimal first decision TI[t, k+1]
                            @debug "Adding penalization..."
                            upd_penalization!(subproblems[b+1], TIbar, kappas[i])
                            # @debug "Setting start values..." # disabled for debug
                            # set_start_value.(model(subproblems[b+1])[:TI], TIbar)
                            @debug "Optimizing..."
                            @debug begin
                                @debug "TR" problem[:TR]
                                show(model(subproblems[b+1]))
                                print(model(subproblems[b+1]))
                                write_to_file(model(subproblems[b+1]), "model.lp")
                            end
                            # @debug begin
                            #     # error("Debug stop")
                            #     undo = relax_integrality(model(subproblems[b+1]))
                            #     # set_optimizer(model(subproblems[b+1]), CDD.Optimizer)
                            #     cons = find_problematic_constraint!(model(subproblems[b+1]))
                            #     open("tmp_infeasibilities.txt", "w") do io
                            #         show(io, "text/plain", cons)
                            #     end
                            #     error("Stop")
                            # end
                            optimize!(model(subproblems[b+1]))
                            error("Stop")
                            # upd_valueTI!(subproblems[t])
                            TIvalues[i, :, :] .= value.(model(subproblems[b+1])[:TI])
                            optsolpertruck[i][:S] = copy(value.(model(subproblems[b+1])[:S]))
                            # optsolpertruck[i][:SG] = copy(value.(model(subproblems[b+1])[:SG]))
                            # optsolpertruck[i][:SL] = copy(value.(model(subproblems[b+1])[:SL]))
                            # optsolpertruck[i][:SW] = copy(value.(model(subproblems[b+1])[:SW]))
                            # optsolpertruck[i][:SH] = copy(value.(model(subproblems[b+1])[:S])) * problem[:IH]
                            optsolpertruck[i][:SXo] = copy(value.(model(subproblems[b+1])[:SXo]))
                            optsolpertruck[i][:SXe] = copy(value.(model(subproblems[b+1])[:SXe]))
                            optsolpertruck[i][:SYo] = copy(value.(model(subproblems[b+1])[:SYo]))
                            optsolpertruck[i][:SYe] = copy(value.(model(subproblems[b+1])[:SYe]))
                            optsolpertruck[i][:SZe] = copy(value.(model(subproblems[b+1])[:SZe]))
                            optsolpertruck[i][:SO] = copy(value.(model(subproblems[b+1])[:SO]))
                            optsolpertruck[i][:SU] = copy(value.(model(subproblems[b+1])[:SU]))
                            optsolpertruck[i][:SK] = copy(value.(model(subproblems[b+1])[:SK]))
                            optsolpertruck[i][:ST] = i
                        end
                    end
                end
                pt += batchsize
            end
            # Update the mean first decisions:
            # TIbar[k+1] = sum([pi[t] * TI[t, k+1] for t in bold_T])
            # TIbar = sum([value.(model(subproblems[t])[:TI]) for t in 1:nbchosentrucks])
            TIbar = sum([TIvalues[t, :, :] for t in 1:nbchosentrucks])

            # Update the multipliers by
            for t in 1:nbchosentrucks
                kappas[t] = kappas[t] .+ delta * sum(TIvalues[t, :, :] .- TIbar)
            end
            TIequality = are_all_TI_equal(TIvalues, TIbar, nbchosentrucks, eps)
        end
    end


    return TIbar, optsolpertruck
end

function TSIProblem(optimizer, instancepath::String)
    item_productcodes, truckdict, supplierdict, supplierdockdict, plantdict, 
    plantdockdict, nbplannedtrucks, nbitems, nbsuppliers, nbsupplierdocks, nbplants, 
    nbplantdocks, TE_P, TL_P, TW_P, TH_P, TKE_P, TGE_P, TDA_P, TDE_P, TU_P, TP_P, TK_P, 
    TG_P, TR_P, IU, IP, IK, IPD, IS, _IO, IL, IW, IH, IDL, IDE, stackabilitycodedict, nbtrucks, TE, 
    TL, TW, TH, TKE, TGE, TDA, TDE, TU, TP, TK, TG, TR, TID, reverse_truckdict, truckindices,
    costinventory, costtransportation, costextratruck, timelimit = loadinstance(instancepath)
    

    return TSIProblem(optimizer, 
    item_productcodes, truckdict, supplierdict, supplierdockdict, plantdict, 
    plantdockdict, nbplannedtrucks, nbitems, nbsuppliers, nbsupplierdocks, nbplants, 
    nbplantdocks, TE_P, TL_P, TW_P, TH_P, TKE_P, TGE_P, TDA_P, TDE_P, TU_P, TP_P, TK_P, 
    TG_P, TR_P, IU, IP, IK, IPD, IS, _IO, IL, IW, IH, IDL, IDE, stackabilitycodedict, nbtrucks, TE, 
    TL, TW, TH, TKE, TGE, TDA, TDE, TU, TP, TK, TG, TR, TID, reverse_truckdict, truckindices,
    costinventory, costtransportation, costextratruck, timelimit
    )

end

getmodel(pb::TSIProblem) = pb.model

function TSIProblem(optimizer, 
    item_productcodes, truckdict, supplierdict, supplierdockdict, plantdict, 
    plantdockdict, nbplannedtrucks, nbitems, nbsuppliers, nbsupplierdocks, nbplants, 
    nbplantdocks, TE_P, TL_P, TW_P, TH_P, TKE_P, TGE_P, TDA_P, TDE_P, TU_P, TP_P, TK_P, 
    TG_P, TR_P, IU, IP, IK, IPD, IS, _IO, IL, IW, IH, IDL, IDE, stackabilitycodedict, nbtrucks, TE, 
    TL, TW, TH, TKE, TGE, TDA, TDE, TU, TP, TK, TG, TR, TID, reverse_truckdict, truckindices,
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
    pb[:TDE_P] = TDE_P
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
    pb[:TDE] = TDE
    pb[:TU] = TU
    pb[:TP] = TP
    pb[:TK] = TK
    pb[:TG] = TG
    pb[:TR] = TR
    pb[:TID] = TID
    pb[:reverse_truckdict] = reverse_truckdict
    pb[:truckindices] = truckindices
    return pb

end


function changetruck!(t, subproblem::Subproblem, chosentrucks; changeS=false)
    
    ## Create model
    submodel = model(subproblem)
    nbcandidateitems = sum(subproblem[:TR][t, :])
    # nbcandidateitems = max([sum(problem[:TR][t2, :]) for t2 in 1:nbtrucks]...)
    nbstacks = nbcandidateitems
    nbitems = subproblem[:nbitems]
    nbtrucks = subproblem[:nbtrucks]
    nbchosentrucks = length(chosentrucks)
    nbplannedtrucks = subproblem[:nbplannedtrucks]
    # nbplants = subproblem[:nbplants]
    nbplantdocks = subproblem[:nbplantdocks]
    nbsuppliers = subproblem[:nbsuppliers]
    nbsupplierdocks = subproblem[:nbsupplierdocks]
    # nbextratrucks = nbtrucks - nbplannedtrucks
    nbextratrucks = length(filter(x -> x in 1:nbplannedtrucks, chosentrucks))
    subproblem.t = t
    # subproblem = Subproblem(t, subproblem, submodel, optimizer, nbstacks, Matrix{Union{Missing, Bool}}(missing, nbtrucks, nbtrucks))

    if changeS
        @info "Replacing S..."
        delete.(submodel, submodel[:S])
        JuMP.unregister.(submodel, :S)
        @variable(submodel, S[1:nbstacks, 1:nbcandidateitems], lower_bound = 0, upper_bound = 1, Bin)
    end
    @info "Computing parameters..."
    MI4 = [[max(subproblem[:TU][:, j]...) for j in 1:first(size(subproblem[:TU][1, :]))] for j in 1:first(size(submodel[:TI][:, 1]))]
    MI4 = map(Int, (max(subproblem[:TU][:, i]...) for i = 1:first(size(subproblem[:TU][1, :])), j = 1:first(size(submodel[:TI][:, 1]))))

    MZ = max(subproblem[:IS]...) + 1.0

    MQ = max(subproblem[:IU]...) + 1.0

    # MH = max(IP...) + 1.0

    MV = max(subproblem[:IK]...) + 1.0

    MW = max(subproblem[:IPD]...) + 1.0

    # MG = max(skipmissing(subproblem[:_IO])...) + 1.0
    MG = 5.0

    MDL =  max(subproblem[:IL]...) + 1.0
    MDW =  max(subproblem[:IW]...) + 1.0

    # Meta = nbtrucks * Mtau/10

    MTL = Matrix{Float64}(undef, nbchosentrucks, 1)
    MTW = Matrix{Float64}(undef, nbchosentrucks, 1)

    MTL = max(subproblem[:TL]...) + 1.0
    MTW = max(subproblem[:TW]...) + 1.0
    # if true
    # # if typeof(MTL) != Vector{T} where {T <: Real}
    #     throw(TypeError(MTL, "MTL must be of type Vector{T} where {T <: Real}", Vector{T} where {T <: Real}, typeof(MTL)))
    # end

    # MTE = max(skipmissing(subproblem[:TE])...)

    MTKE = max(subproblem[:TKE]...)

    MTGE = max(subproblem[:TGE]...)

    MTW = max(subproblem[:TW]...) + 1.0

    # MPsi = nbitems

    Mlambda = 2.0 * subproblem[:TL][t] + 1.0
    Mlambda = 2.0 * max(subproblem[:TL]...) + 1.0

    Mzeta = nbitems

    # SZo = 0.0

    # MST = 2.0

    # MOmega = 2.0

    # Mmu = 2.0

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
        delete.(submodel, submodel[:cZetaT2])
        JuMP.unregister.(submodel, :cZetaT2)
        @constraint(submodel, cZetaT2, submodel[:zetaT] * Mzeta .>= submodel[:TI][1:nbplannedtrucks, :] * vones(Int8, nbitems))
        delete.(submodel, submodel[:cZetaT3])
        JuMP.unregister.(submodel, :cZetaT3)
        @constraint(submodel, cZetaT3, -(1 .- submodel[:zetaT]) * Mzeta + 1 .<= submodel[:TI][1:nbplannedtrucks, :] * vones(Int8, nbitems))
else
        @info "Replacing cZetaE2..."
        delete.(submodel, submodel[:cZetaE2])
        JuMP.unregister.(submodel, :cZetaE2)
        @constraint(submodel, cZetaE2, submodel[:zetaE] * Mzeta .>= submodel[:TI][nbplannedtrucks+1:end, :] * vones(Int8, nbitems))
        delete.(submodel, submodel[:cZetaE3])
        JuMP.unregister.(submodel, :cZetaE3)
        @constraint(submodel, cZetaE3, -(1 .- submodel[:zetaE]) * Mzeta + 1 .<= submodel[:TI][nbplannedtrucks+1:end, :] * vones(Int8, nbitems))
end
    @info "Replacing cTI_TR..."
    delete.(submodel, submodel[:cTI_TR])
    JuMP.unregister.(submodel, :cTI_TR)
    @constraint(submodel, cTI_TR, submodel[:TI][t, :] .<= subproblem[:TR][t, :])

    icandidates = findall((x) -> x == 1, subproblem[:TR][t, :])
    @info "Replacing cTI_1_1..."
    delete.(submodel, submodel[:cTI_1_1])
    JuMP.unregister.(submodel, :cTI_1_1)
    @constraint(submodel, cTI_1_1, transpose(submodel[:TI])[filter(x -> x in icandidates, 1:nbitems), :] * vones(Int8, nbchosentrucks) + submodel[:GI][filter(x -> x in icandidates, 1:nbitems)] .<= vones(Int8, size(transpose(submodel[:TI])[filter(x -> x in icandidates, 1:nbitems), :], 1)))
    # @debug "cTI_1_1" cTI_1_1[icandidates[1], :]
    @info "Replacing cS_TI..."
    delete.(submodel, submodel[:cS_TI])
    JuMP.unregister.(submodel, :cS_TI)
    @constraint(submodel, cS_TI, transpose(submodel[:S]) * vones(Int8, nbstacks) .== submodel[:TI][t, filter(x -> x in icandidates, 1:nbitems)])

    if changeS
        @info "Replacing cZ_S_MZ..."
        delete.(submodel, submodel[:cZ_S_MZ])
        JuMP.unregister.(submodel, :cZ_S_MZ)
        @constraint(submodel, cZ_S_MZ, submodel[:Z] .<= submodel[:S] * MZ)
    end
    # @info "Replacing cS_TI_MS..."
    # delete.(submodel, submodel[:cS_TI_MSleft])
    # JuMP.unregister.(submodel, :cS_TI_MSleft)
    # @constraint(submodel, cS_TI_MSleft, -(vones(Int8, nbcandidateitems) - submodel[:TI][t, filter(x -> x in icandidates, 1:nbitems)]) * MS .<= transpose(submodel[:S]) * vones(Int8, nbcandidateitems) - vones(Int8, nbcandidateitems))
    # delete.(submodel, submodel[:cS_TI_MSright])
    # JuMP.unregister.(submodel, :cS_TI_MSright)
    # @constraint(submodel, cS_TI_MSright, transpose(submodel[:S]) * vones(Int8, nbcandidateitems) - vones(Int8, nbcandidateitems) .<= (vones(Int8, nbcandidateitems) - submodel[:TI][t, filter(x -> x in icandidates, 1:nbitems)]) * MS)

    @info "Replacing cS_IS_Z..."
    delete.(submodel, submodel[:cS_IS_Z])
    JuMP.unregister.(submodel, :cS_IS_Z)
    @constraint(submodel, cS_IS_Z, submodel[:S] * subproblem[:IS][filter(x -> x in icandidates, 1:nbitems)] .== submodel[:Z] * vones(Int8, size(submodel[:Z])[1]))

    if changeS
        @info "Replacing cZ_SS_S..."
        delete.(submodel, submodel[:cZ_SS_Sleft])
        JuMP.unregister.(submodel, :cZ_SS_Sleft)
        @constraint(submodel, cZ_SS_Sleft, -MZ * (-submodel[:S] .+ 1) .<= submodel[:Z] .- hcat([submodel[:SS] for i in 1:size(submodel[:Z])[2]]...))
        delete.(submodel, submodel[:cZ_SS_Sright])
        JuMP.unregister.(submodel, :cZ_SS_Sright)
        @constraint(submodel, cZ_SS_Sright, submodel[:Z] .- hcat([submodel[:SS] for i in 1:size(submodel[:Z])[2]]...) .<= MZ * (-submodel[:S] .+ 1))
    end

    @info "Replacing cS_IU_Q..."
    delete.(submodel, submodel[:cS_IU_Q])
    JuMP.unregister.(submodel, :cS_IU_Q)
    @constraint(submodel, cS_IU_Q, submodel[:S] * subproblem[:IU][filter(x -> x in icandidates, 1:nbitems)] .== submodel[:Q])
    
    if changeS
        @info "Replacing cQ_SU_S..."
        delete.(submodel, submodel[:cQ_SU_Sleft])
        JuMP.unregister.(submodel, :cQ_SU_Sleft)
        @constraint(submodel, cQ_SU_Sleft, -MQ * (1 .- submodel[:SU]) .<= submodel[:Q] .- hcat([submodel[:S] * vones(Int8, size(submodel[:S], 2)) for i in 1:size(submodel[:Q], 2)]...))
        delete.(submodel, submodel[:cQ_SU_Sright])
        JuMP.unregister.(submodel, :cQ_SU_Sright)
        @constraint(submodel, cQ_SU_Sright, submodel[:Q] .- hcat([submodel[:S] * vones(Int8, size(submodel[:S], 2)) for i in 1:size(submodel[:Q], 2)]...) .<= MQ * (1 .- submodel[:SU]))
    end

    @info "Replacing cS_IK_V..."
    delete.(submodel, submodel[:cS_IK_V])
    JuMP.unregister.(submodel, :cS_IK_V)
    @constraint(submodel, cS_IK_V, submodel[:S] * subproblem[:IK][filter(x -> x in icandidates, 1:nbitems)] .== submodel[:V])

    if changeS
        @info "Adding cV_SK_S..."
        delete.(submodel, submodel[:cV_SK_Sleft])
        JuMP.unregister.(submodel, :cV_SK_Sleft)
        @constraint(submodel, cV_SK_Sleft, -MV * (-submodel[:SK] .+ 1) .<= submodel[:V] .- hcat([submodel[:S] * vones(Int8, size(submodel[:S], 2)) for i in 1:size(submodel[:V], 2)]...))
        delete.(submodel, submodel[:cV_SK_Sright])
        JuMP.unregister.(submodel, :cV_SK_Sright)
        @constraint(submodel, cV_SK_Sright, submodel[:V] .- hcat([submodel[:S] * vones(Int8, size(submodel[:S], 2)) for i in 1:size(submodel[:V], 2)]...) .<= MV * (-submodel[:SK] .+ 1))
    end

    @info "Replacing cS_IPD_W..."
    delete.(submodel, submodel[:cS_IPD_W])
    JuMP.unregister.(submodel, :cS_IPD_W)
    @constraint(submodel, cS_IPD_W, submodel[:S] * subproblem[:IPD][filter(x -> x in icandidates, 1:nbitems)] .== submodel[:W])

    if changeS
        @info "Replacing cW_SG_S..."
        delete.(submodel, submodel[:cW_SPD_Sleft])
        JuMP.unregister.(submodel, :cW_SPD_Sleft)
        @constraint(submodel, cW_SPD_Sleft, -MW * (-submodel[:SG] .+ 1) .<= submodel[:W] .- hcat([submodel[:S] * vones(Int8, size(submodel[:S], 2)) for i in 1:size(submodel[:W], 2)]...))
        delete.(submodel, submodel[:cW_SPD_Sright])
        JuMP.unregister.(submodel, :cW_SPD_Sright)
        @constraint(submodel, cW_SPD_Sright, submodel[:W] .- hcat([submodel[:S] * vones(Int8, size(submodel[:S], 2)) for i in 1:size(submodel[:W], 2)]...) .<= MW * (-submodel[:SG] .+ 1))

        @info "Replacing cGl_S..."
        delete.(submodel, submodel[:cGl_S])
        JuMP.unregister.(submodel, :cGl_S)
        @constraint(submodel, cGl_S, submodel[:Gl] .<= submodel[:S] * MG)
        @info "Replacing cGr_S..."
        delete.(submodel, submodel[:cGr_S])
        JuMP.unregister.(submodel, :cGr_S)
        @constraint(submodel, cGr_S, submodel[:Gr] .<= submodel[:S] * MG)

        @info "Replacing cGr_SO_S..."
        delete.(submodel, submodel[:cGr_SO_Sleft])
        JuMP.unregister.(submodel, :cGr_SO_Sleft)
        @constraint(submodel, cGr_SO_Sleft, -MG * (-submodel[:S] .+ 1) .<= submodel[:Gr] .- hcat([submodel[:SO] for i in 1:nbcandidateitems]...))
        delete.(submodel, submodel[:cGr_SO_Sright])
        JuMP.unregister.(submodel, :cGr_SO_Sright)
        @constraint(submodel, cGr_SO_Sright, submodel[:Gr] .- hcat([submodel[:SO] for i in 1:nbcandidateitems]...) .<= MG * (-submodel[:S] .+ 1))
    end

    @info "Replacing cGl_IOV_S..."
    delete.(submodel, submodel[:cGl_IOV_Sleft])
    JuMP.unregister.(submodel, :cGl_IOV_Sleft)
    @constraint(submodel, cGl_IOV_Sleft, -MG * (-submodel[:S] .+ 1) .<= submodel[:Gl] .- vcat([transpose(submodel[:IOV][filter(x -> x in icandidates, 1:nbitems)]) for i in 1:nbstacks]...))
    delete.(submodel, submodel[:cGl_IOV_Sright])
    JuMP.unregister.(submodel, :cGl_IOV_Sright)
    @constraint(submodel, cGl_IOV_Sright, submodel[:Gl] .- vcat([transpose(submodel[:IOV][filter(x -> x in icandidates, 1:nbitems)]) for i in 1:nbstacks]...) .<= MG * (-submodel[:S] .+ 1))

    if changeS
        @info "Replacing cDL_S..."
        delete.(submodel, submodel[:cDL_S])
        JuMP.unregister.(submodel, :cDL_S)
        @constraint(submodel, cDL_S, submodel[:DL] .<= submodel[:S] * MDL)
        @info "Replacing cDL_SL..."
        delete.(submodel, submodel[:cDL_SLleft])
        JuMP.unregister.(submodel, :cDL_SLleft)
        @constraint(submodel, cDL_SLleft, -MDL * (-submodel[:S] .+ 1) .<= submodel[:DL] - hcat([submodel[:SL] for i in 1:size(submodel[:DL])[2]]...))
        delete.(submodel, submodel[:cDL_SLright])
        JuMP.unregister.(submodel, :cDL_SLright)
        @constraint(submodel, cDL_SLright, submodel[:DL] - hcat([submodel[:SL] for i in 1:size(submodel[:DL])[2]]...) .<= MDL * (-submodel[:S] .+ 1))
    end
    @info "Replacing cDL_S_IL..."
    delete.(submodel, submodel[:cDL_S_IL])
    JuMP.unregister.(submodel, :cDL_S_IL)
    @constraint(submodel, cDL_S_IL, submodel[:DL] * vones(Int8, size(submodel[:S], 1)) .== submodel[:S] * subproblem[:IL][filter(x -> x in icandidates, 1:nbitems)])
    
    if changeS
        @info "Replacing cDW_S..."
        delete.(submodel, submodel[:cDW_S])
        JuMP.unregister.(submodel, :cDW_S)
        @constraint(submodel, cDW_S, submodel[:DW] .<= submodel[:S] * MDW)
        @info "Replacing cDW_SW..."
        delete.(submodel, submodel[:cDW_SWleft])
        JuMP.unregister.(submodel, :cDW_SWleft)
        @constraint(submodel, cDW_SWleft, -MDW * (-submodel[:S] .+ 1) .<= submodel[:DW] - hcat([submodel[:SW] for i in 1:size(submodel[:DW])[2]]...))
        delete.(submodel, submodel[:cDW_SWright])
        JuMP.unregister.(submodel, :cDW_SWright)
        @constraint(submodel, cDW_SWright, submodel[:DW] - hcat([submodel[:SW] for i in 1:size(submodel[:DW])[2]]...) .<= MDW * (-submodel[:S] .+ 1))
    end


    @info "Replacing cDW_S_IW..."
    delete.(submodel, submodel[:cDW_S_IW])
    JuMP.unregister.(submodel, :cDW_S_IW)
    @constraint(submodel, cDW_S_IW, submodel[:DW] * vones(Int8, size(submodel[:S], 1)) .== submodel[:S] * subproblem[:IW][filter(x -> x in icandidates, 1:nbitems)])

    @info "Replacing cSZe_S_IH..."
    delete.(submodel, submodel[:cSZe_S_IH])
    JuMP.unregister.(submodel, :cSZe_S_IH)
    @constraint(submodel, cSZe_S_IH, submodel[:SZe] .== submodel[:S]* subproblem[:IH][filter(x -> x in icandidates, 1:nbitems)])

    @info "Replacing cSXe_ST_TL..."
    delete.(submodel, submodel[:cSXe_ST_TL])
    JuMP.unregister.(submodel, :cSXe_ST_TL)
    @constraint(submodel, cSXe_ST_TL, submodel[:SXe] .<= subproblem[:TL][t])
    @info "Replacing cSYe_ST_TW..."
    delete.(submodel, submodel[:cSYe_ST_TW])
    JuMP.unregister.(submodel, :cSYe_ST_TW)
    @constraint(submodel, cSYe_ST_TW, submodel[:SYe] .<= subproblem[:TW][t])
    @info "Replacing cSZe_ST_TH..."
    delete.(submodel, submodel[:cSZe_ST_TH])
    JuMP.unregister.(submodel, :cSZe_ST_TH)
    @constraint(submodel, cSZe_ST_TH, submodel[:SZe] .<= subproblem[:TH][t])

    @info "Replacing cXi1SU_Xi2SU..."
    notmissingTE = filter(x -> !ismissing(subproblem[:TE][t, x]), 1:nbsuppliers)
    delete.(submodel, submodel[:cXi1SU_Xi2SU])
    JuMP.unregister.(submodel, :cXi1SU_Xi2SU)
    @constraint(submodel, cXi1SU_Xi2SU, Xi1 * submodel[:SU][:, notmissingTE] * subproblem[:TE][t, notmissingTE] .<= Xi2 * submodel[:SU][:, notmissingTE] * subproblem[:TE][t, notmissingTE])

    @debug begin
        @debug "cXi1SU_Xi2SU" cXi1SU_Xi2SU
        sleep(20)
    end

    @info "Replacing cXi2SK_Xi1SK..."
    delete.(submodel, submodel[:cXi2SK_Xi1SK])
    JuMP.unregister.(submodel, :cXi2SK_Xi1SK)
    notmissingTKE = filter(x -> !ismissing(subproblem[:TKE][t, x]), 1:nbsupplierdocks)
    @constraint(submodel, cXi2SK_Xi1SK, Xi2*submodel[:SK][:, notmissingTKE]*subproblem[:TKE][t, notmissingTKE] .>= Xi1*submodel[:SK][:, notmissingTKE]*subproblem[:TKE][t, notmissingTKE] - (-submodel[:r] .+ 1) * MTKE)

    @info "Replacing cXi1SK_Xi2SK_chi..."
    delete.(submodel, submodel[:cXi1SK_Xi2SK_chi])
    JuMP.unregister.(submodel, :cXi1SK_Xi2SK_chi)
    @constraint(submodel, cXi1SK_Xi2SK_chi, Xi1*submodel[:SK][:, notmissingTKE]*subproblem[:TKE][t, notmissingTKE] - Xi2*submodel[:SK][:, notmissingTKE]*subproblem[:TKE][t, notmissingTKE] .>= submodel[:chi]*epsilon - (-submodel[:sigma2] .+ 1)*MTKE)
    @info "Replacing cXi2SK_Xi1SK_chi..."
    delete.(submodel, submodel[:cXi2SK_Xi1SK_chi])
    JuMP.unregister.(submodel, :cXi2SK_Xi1SK_chi)
    @constraint(submodel, cXi2SK_Xi1SK_chi, Xi2*submodel[:SK][:, notmissingTKE]*subproblem[:TKE][t, notmissingTKE] - Xi1*submodel[:SK][:, notmissingTKE]*subproblem[:TKE][t, notmissingTKE] .>= (-submodel[:chi] .+ 1)*epsilon - (-submodel[:sigma2] .+ 1)*MTKE)

    @info "Replacing cXi2SG_Xi1SG..."
    delete.(submodel, submodel[:cXi2SG_Xi1SG])
    JuMP.unregister.(submodel, :cXi2SG_Xi1SG)
    notmissingTGE = filter(x -> !ismissing(subproblem[:TGE][t, x]), 1:nbplantdocks)
    @constraint(submodel, cXi2SG_Xi1SG, Xi2*submodel[:SG][:, notmissingTGE]*subproblem[:TGE][t, notmissingTGE] .>= Xi1*submodel[:SG][:, notmissingTGE]*subproblem[:TGE][t, notmissingTGE] - (-submodel[:sigma3] .+ 1) * MTGE)
end

function Subproblem(t, problem, optimizer, chosentrucks)
    
    ## Create model
    submodel = Model(optimizer)
    nbcandidateitems = sum(problem[:TR][t, :])
    nbitems = problem[:nbitems]
    nbtrucks = problem[:nbtrucks]
    nbchosentrucks = length(chosentrucks)
    # nbcandidateitems = max([sum(problem[:TR][t2, :]) for t2 in 1:nbtrucks]...)
    nbstacks = nbcandidateitems
    nbplannedtrucks = problem[:nbplannedtrucks]
    nbplants = problem[:nbplants]
    nbplantdocks = problem[:nbplantdocks]
    nbsuppliers = problem[:nbsuppliers]
    nbsupplierdocks = problem[:nbsupplierdocks]
    nbextratrucks = length(filter(x -> x in 1:nbplannedtrucks, chosentrucks))
    subproblem = Subproblem(t, problem, submodel, optimizer, nbstacks, Matrix{Union{Missing, Bool}}(missing, nbtrucks, nbtrucks))

    ## Add variables
    @info "Creating variables..."
    @info "Adding zetaT..."
    @variable(submodel, zetaT[1:nbplannedtrucks] >= 0) 
    @info "Adding zetaE..."
    @variable(submodel, zetaE[1:nbextratrucks] >= 0)

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
    @variable(submodel, TI[1:nbchosentrucks, 1:nbitems], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding GI..."
    @variable(submodel, GI[1:nbitems], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding R..."
    @variable(submodel, R[1:nbitems, 1:nbsuppliers], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding Theta..."
    @variable(submodel, Theta[1:nbitems, 1:nbsuppliers], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding S..."
    @variable(submodel, S[1:nbstacks, 1:nbcandidateitems], lower_bound = 0, upper_bound = 1, Bin)
    @info "Adding Z..."
    @variable(submodel, Z[1:nbstacks, 1:nbcandidateitems], lower_bound = 0, upper_bound = 1, Bin)

    @debug "nbstacks" nbstacks
    @debug "nbtrucks" nbtrucks
    @debug "nbchosentrucks" nbchosentrucks
    @debug "nbitems" nbitems
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

    MG = 5.0

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

    Mzeta = nbitems

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


    ### DEBUG begin
    ## Add constraints
    @info "Adding constraints..."
    if t <= nbplannedtrucks
        @info "Adding cZetaT2..."
        # with_logger(ConsoleLogger(stdout, Logging.Debug, show_limited=true)) do
        #     @debug "-reshape(TI[t, :], nbitems, 1)" -reshape(TI[t, :], nbitems, 1)
        #     @debug "vones(Int8, nbitems)" vones(Int8, nbitems)
        # end
        @constraint(submodel, cZetaT2, zetaT * Mzeta .>= TI[1:nbplannedtrucks, :] * vones(Int8, nbitems))
        @info "Adding cZetaT3..."
        @constraint(submodel, cZetaT3, -(1 .- zetaT) * Mzeta .+ 1 .<= TI[1:nbplannedtrucks, :] * vones(Int8, nbitems))
    else
        @info "Adding cZetaE2..."
        @constraint(submodel, cZetaE2, zetaE * Mzeta .>= TI[nbplannedtrucks+1:end, :] * vones(Int8, nbitems))
        @info "Adding cZetaE3..."
        @constraint(submodel, cZetaE3, -(1 .- zetaE) * Mzeta .+ 1 .<= TI[nbplannedtrucks+1:end, :] * vones(Int8, nbitems))
    end
    @info "Adding cTI_TR..."
    @constraint(submodel, cTI_TR, TI[t, :] .<= problem[:TR][t, :])

    icandidates = findall((x) -> x == 1, problem[:TR][t, :])
    # @debug icandidates
    @info "Adding cTI_1_1..."
    # no more than one truck per candidate item
    # with_logger(ConsoleLogger(stdout, Logging.Debug, show_limited=true)) do
    #     @debug "transpose(TI)[filter(x -> x in icandidates, 1:nbitems), :]" transpose(TI)[filter(x -> x in icandidates, 1:nbitems), :]
    #     @debug "vones(Int8, nbtrucks)" vones(Int8, nbtrucks)
    # end
    @constraint(submodel, cTI_1_1, transpose(TI)[filter(x -> x in icandidates, 1:nbitems), :] * vones(Int8, nbchosentrucks) + GI[filter(x -> x in icandidates, 1:nbitems)] .<= vones(Int8, size(transpose(TI)[filter(x -> x in icandidates, 1:nbitems), :], 1)))
    # @debug "cTI_1_1" cTI_1_1[icandidates[1], :]
    @info "Adding cS_TI..."
    # with_logger(ConsoleLogger(stdout, Logging.Debug, show_limited=true)) do 
    #     @debug "S" S
    #     @debug "vones(Int8, length(icandidates))" vones(Int8, length(icandidates))
    #     @debug "S * vones(Int8, length(icandidates))" S * vones(Int8, length(icandidates))
    # end
    @constraint(submodel, cS_TI, transpose(S) * vones(Int8, nbstacks) .== TI[t, filter(x -> x in icandidates, 1:nbitems)])
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

    # @info "Adding cS_TI_MS..."
    # @constraint(submodel, cS_TI_MSleft, -(vones(Int8, nbcandidateitems) - TI[t, filter(x -> x in icandidates, 1:nbitems)]) * MS .<= transpose(S) * vones(Int8, nbcandidateitems) - vones(Int8, nbcandidateitems))
    # @constraint(submodel, cS_TI_MSright, transpose(S) * vones(Int8, nbcandidateitems) - vones(Int8, nbcandidateitems) .<= (vones(Int8, nbcandidateitems) - TI[t, filter(x -> x in icandidates, 1:nbitems)]) * MS)

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

    @info "Adding cGl_Gr..."
    # @constraint(submodel, cGl_Gr, Gl * vones(Int8, size(Gl)[1]) .== Gr * vones(Int8, size(Gr)[1])) # I don't understand this
    @constraint(submodel, cGl_Gr, Gl .== Gr)

    @info "Adding cGr_SO_S..."
    @debug begin
        @debug "SO " SO 
        # sleep(10)
        @debug "hcat([SO for i in 1:size(Gr)[2]]...)" hcat([SO for i in 1:nbcandidateitems]...)
        # sleep(10)
        @debug "Gr" Gr
        # sleep(30)
    end
    @constraint(submodel, cGr_SO_Sleft, -MG * (-S .+ 1) .<= Gr .- hcat([SO for i in 1:nbcandidateitems]...))
    @constraint(submodel, cGr_SO_Sright, Gr .- hcat([SO for i in 1:nbcandidateitems]...) .<= MG * (-S .+ 1))
    @info "Adding cGl_IOV_S..."
    @debug begin
        @debug "IOV " IOV 
        # sleep(10)
        @debug "vcat([transpose(IOV[filter(x -> x in icandidates, 1:nbitems)]) for i in 1:size(Gl)[2]]...)" vcat([transpose(IOV[filter(x -> x in icandidates, 1:nbitems)]) for i in 1:size(Gl)[2]]...)
        # sleep(10)
        @debug "Gl" Gl
        # sleep(30)
    end
    @constraint(submodel, cGl_IOV_Sleft, -MG * (-S .+ 1) .<= Gl .- vcat([transpose(IOV[filter(x -> x in icandidates, 1:nbitems)]) for i in 1:nbstacks]...))
    @constraint(submodel, cGl_IOV_Sright, Gl .- vcat([transpose(IOV[filter(x -> x in icandidates, 1:nbitems)]) for i in 1:nbstacks]...) .<= MG * (-S .+ 1))

    @info "Adding cDL_S..."
    @constraint(submodel, cDL_S, DL .<= S * MDL)
    ## debug
    @info "Adding cDL_SL..."
    @constraint(submodel, cDL_SLleft, -MDL * (-S .+ 1) .<= DL - hcat([SL for i in 1:size(DL)[2]]...))
    @constraint(submodel, cDL_SLright, DL - hcat([SL for i in 1:size(DL)[2]]...) .<= MDL * (-S .+ 1))
    @info "Adding cDL_S_IL..."
    @constraint(submodel, cDL_S_IL, DL * vones(Int8, size(S, 1)) .== S * problem[:IL][filter(x -> x in icandidates, 1:nbitems)])
    @info "Adding cDW_S..."
    @constraint(submodel, cDW_S, DW .<= S * MDW)
    ## debug
    @info "Adding cDW_SW..."
    @constraint(submodel, cDW_SWleft, -MDW * (-S .+ 1) .<= DW - hcat([SW for i in 1:size(DW)[2]]...))
    @constraint(submodel, cDW_SWright, DW - hcat([SW for i in 1:size(DW)[2]]...) .<= MDW * (-S .+ 1))
    @info "Adding cDW_S_IW..."
    @constraint(submodel, cDW_S_IW, DW * vones(Int8, size(S, 1)) .== S * problem[:IW][filter(x -> x in icandidates, 1:nbitems)])

    ## debug
    @info "Adding cSXe_SXo_SL_SO..."
    # @constraint(submodel, cSXe_SXo_SL_SO, SXe - SXo .== SL + SO * MTL)
    @constraint(submodel, cSXe_SXo_SL_SOleft, SXe - SXo - SL .<= SO * MTL)
    @constraint(submodel, cSXe_SXo_SL_SOright, SXe - SXo - SL .>= -SO * MTL)
    @info "Adding cSYe_SYo_SW_SO..."
    # @constraint(submodel, cSYe_SYo_SW_SO, SYe - SYo .== SW + SO * MTW)
    @constraint(submodel, cSYe_SYo_SW_SOleft, SYe - SYo - SW .<= SO * MTW)
    @constraint(submodel, cSYe_SYo_SW_SOright, SYe - SYo - SW .>= -SO * MTW)
    @info "Adding cSXe_SXo_SW_SO..."
    # @constraint(submodel, cSXe_SXo_SW_SO, SXe - SXo .== SW + (-SO .+ 1) * MTW)
    @constraint(submodel, cSXe_SXo_SW_SOleft, SXe - SXo - SW .<= (-SO .+ 1) * MTW)
    @constraint(submodel, cSXe_SXo_SW_SOright, SXe - SXo - SW .>= -(-SO .+ 1) * MTW)
    @info "Adding cSYe_SYo_SL_SO..."
    # @constraint(submodel, cSYe_SYo_SL_SO, SYe - SYo .== SL + (-SO .+ 1) * MTL)
    @constraint(submodel, cSYe_SYo_SL_SOleft, SYe - SYo - SL .<= (-SO .+ 1) * MTL)
    @constraint(submodel, cSYe_SYo_SL_SOright, SYe - SYo - SL .>= -(-SO .+ 1) * MTL)
    @info "Adding cSZe_S_IH..."
    @constraint(submodel, cSZe_S_IH, SZe .== S* problem[:IH][filter(x -> x in icandidates, 1:nbitems)])
    @info "Adding cSXe_ST_TL..."
    @constraint(submodel, cSXe_ST_TL, SXe .<= problem[:TL][t])
    @info "Adding cSYe_ST_TW..."
    @constraint(submodel, cSYe_ST_TW, SYe .<= problem[:TW][t])
    @info "Adding cSZe_ST_TH..."
    @constraint(submodel, cSZe_ST_TH, SZe .<= problem[:TH][t])
    #####

    @info "Adding cSXo_SXo..."
    # with_logger(ConsoleLogger(stdout, Logging.Debug, show_limited=true)) do
    #     @debug  "SXo[1:end-1]" SXo[1:end-1]
    #     @debug "SXo[2:end]" SXo[2:end]
    # end
    # @constraint(submodel, cSXo_SXo, (vcat(hcat([1], falses(1, size(SXo)[1]-1)), I(size(SXo)[1])) * SXo)[2:end] .<= SXo)
    
    ## debug
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

    ## debug
    @constraint(submodel, cXi2SXo_Xi1SXe_betaM_betaP, (Xi2 * SXo) - (Xi1 * SXe) - betaM + betaP .== -epsilon * vones(Float64, size(Xi1, 1)))

    @info "Adding cbetaM_lambda..."
    # with_logger(ConsoleLogger(stdout, Logging.Debug, show_limited=true)) do
    #     @debug "betaM" betaM
    #     @debug "lambda" lambda
    #     @debug "Mlambda" Mlambda 
    # end

    ## debug
    @constraint(submodel, cbetaM_lambda, betaM .<= lambda .* Mlambda)

    @info "Adding cbetaP_lambda..."

    ## debug
    @constraint(submodel, cbetaP_lambda, betaP .<= (-lambda .+ 1)*Mlambda)

    @info "Adding cmu_betaM..."

    ## debug
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

    ## debug
    @constraint(submodel, cXi1SYe_Xi2SYo, Xi1 * SYe .<= Xi2 * SYo + xi * MTW + (-mu .+ 1) * MTW)
    @info "Adding cXi2SYe_Xi1SYo..."

    ## debug
    @constraint(submodel, cXi2SYe_Xi1SYo, Xi2 * SYe .<= Xi1 * SYo + (-xi .+ 1) * MTW + (-mu .+ 1) * MTW)

    @info "Adding cXi1SU_Xi2SU..."
    # with_logger(ConsoleLogger(stdout, Logging.Debug, show_limited=true)) do
    #     @debug "Xi1" Xi1
    #     @debug "SU" SU
    #     @debug "problem[:TE]" problem[:TE]
    #     @debug "problem[:TE][t, :]" problem[:TE][t, :]
    # end
    notmissingTE = filter(x -> !ismissing(problem[:TE][t, x]), 1:nbsuppliers)

    ## debug
    @constraint(submodel, cXi1SU_Xi2SU, Xi1 * SU[:, notmissingTE] * problem[:TE][t, notmissingTE] .<= Xi2 * SU[:, notmissingTE] * problem[:TE][t, notmissingTE])
    # @debug begin
    #     @debug "cXi1SU_Xi2SU" cXi1SU_Xi2SU
    #     @debug "notmissingTE" notmissingTE
    #     @debug "Xi1" Xi1
    #     @debug "Xi2" Xi2
    #     @debug "SU[:, notmissingTE]" SU[:, notmissingTE]
    #     @debug "problem[:TE][t, notmissingTE]" problem[:TE][t, notmissingTE]
    #     sleep(20)
    # end
    # error("Debug stop!!")
    @info "Adding cXi1SU_Xi2SU_chi..."

    ## debug
    @constraint(submodel, cXi1SU_Xi2SU_chi, Xi1*SU - Xi2*SU .>= chi * epsilon - r*MTE - (-sigma1 .+ 1) * MTE)

    @info "Adding cXi2SU_Xi1SU_chi..."
    ## debug
    @constraint(submodel, cXi2SU_Xi1SU_chi, Xi2*SU - Xi1*SU .>= (-chi .+ 1) * epsilon - r*MTE - (-sigma1 .+ 1) * MTE)

    @info "Adding cXi2SK_Xi1SK..."
    notmissingTKE = filter(x -> !ismissing(problem[:TKE][t, x]), 1:nbsupplierdocks)

    ## debug
    @constraint(submodel, cXi2SK_Xi1SK, Xi2*SK[:, notmissingTKE]*problem[:TKE][t, notmissingTKE] .>= Xi1*SK[:, notmissingTKE]*problem[:TKE][t, notmissingTKE] - (-r .+ 1) * MTKE)

    @info "Adding cXi1SK_Xi2SK_chi..."
    
    ## debug
    @constraint(submodel, cXi1SK_Xi2SK_chi, Xi1*SK[:, notmissingTKE]*problem[:TKE][t, notmissingTKE] - Xi2*SK[:, notmissingTKE]*problem[:TKE][t, notmissingTKE] .>= chi*epsilon - (-sigma2 .+ 1)*MTKE)
    @info "Adding cXi2SK_Xi1SK_chi..."

    ## debug
    @constraint(submodel, cXi2SK_Xi1SK_chi, Xi2*SK[:, notmissingTKE]*problem[:TKE][t, notmissingTKE] - Xi1*SK[:, notmissingTKE]*problem[:TKE][t, notmissingTKE] .>= (-chi .+ 1)*epsilon - (-sigma2 .+ 1)*MTKE)

    @info "Adding cXi2SG_Xi1SG..."
    notmissingTGE = filter(x -> !ismissing(problem[:TGE][t, x]), 1:nbplantdocks)
    # with_logger(ConsoleLogger(stdout, Logging.Debug, show_limited=true)) do
    #     @debug "SG[:, notmissingTGE]" SG[:, notmissingTGE]
    #     @debug "problem[:TGE][t, notmissingTGE]" problem[:TGE][t, notmissingTGE]
    # end

    ## debug
    @constraint(submodel, cXi2SG_Xi1SG, Xi2*SG[:, notmissingTGE]*problem[:TGE][t, notmissingTGE] .>= Xi1*SG[:, notmissingTGE]*problem[:TGE][t, notmissingTGE] - (-sigma3 .+ 1) * MTGE)

    @info "Adding csigma1_sigma2_sigma3..."

    ## debug
    @constraint(submodel, csigma1_sigma2_sigma3, sigma1 + sigma2 + sigma3 .>= 1)

    #### DEBUG end
    #### DEBUG begin
    # @info "Adding IOV constraints..."
    # for i in 1:nbitems
    #     if !ismissing(problem[:_IO][i])
    #         submodel[Symbol(string("cIOV[", i, "]"))] = @constraint(submodel, IOV[i] == problem[:_IO][i])
    #     end
    # end
    #### DEBUG end
    # @info "Adding objective function..."
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



    # @objective(submodel, Min, 
    # sum(problem[:costtransportation] * zetaT) + 
    # sum(problem[:costextratruck] * zetaE) + 
    # sum(problem[:costinventory] * (problem[:IDL] - transpose(TI) * problem[:TDA])))

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