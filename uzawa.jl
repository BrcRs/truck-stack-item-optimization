using JuMP
using MathOptInterface

using Base.Threads
using Clp # Only for debug
# using CDD # Only debug
# include("instance_loader.jl")
include("matrix_ops.jl")
include("linear_infeasibilities.jl")
include("progress.jl")
include("subproblem.jl")
include("tsiproblem.jl")

function upd_penalization!(subpb::Subproblem, TIbar, kappa)
    submodel = model(subpb)
    t = truck(subpb)

    nbplannedtrucks = subpb[:nbplannedtrucks]
    nbitems = subpb[:nbitems]
    # @objective(submodel, Min, 
    # sum(problem[:costtransportation] * zetaT) + 
    # sum(problem[:costextratruck] * zetaE) + 
    # sum(problem[:costinventory] * (problem[:IDL] - transpose(TI) * problem[:TDA])))


    # error("Stop")

    # Worst case solution value:
    MGI = subpb[:costtransportation] * nbplannedtrucks + (nbitems - nbplannedtrucks) * subpb[:costextratruck] + subpb[:costinventory] * nbitems * (min(subpb[:IDL]...) - max(subpb[:TDE]...))

    if haskey(JuMP.object_dictionary(submodel), :obj1)
        # delete(submodel, submodel[:obj1])
        JuMP.unregister(submodel, :obj1)
    end
    if haskey(JuMP.object_dictionary(submodel), :obj2)
        # delete(submodel, submodel[:obj2])
        JuMP.unregister(submodel, :obj2)
    end
    if haskey(JuMP.object_dictionary(submodel), :obj3)
        # delete(submodel, submodel[:obj3])
        JuMP.unregister(submodel, :obj3)
    end
    if haskey(JuMP.object_dictionary(submodel), :obj4)
        # delete(submodel, submodel[:obj4])
        JuMP.unregister(submodel, :obj4)
    end
    if haskey(JuMP.object_dictionary(submodel), :obj5)
        # delete(submodel, submodel[:obj5])
        JuMP.unregister(submodel, :obj5)
    end

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
function are_all_TI_equal(TIvalues, TIbar, nbtrucks, eps; retgap=false)
    # @info "Verifying equality..."
    gap = 0
    TIequality = true
    for t in 1:nbtrucks
        diff = sum(TIvalues[t, :, :] == TIbar ? 0.0 : abs.(TIvalues[t, :, :] .- TIbar))
        if diff > eps
            TIequality = false
            # @info "All not equal"
            if !retgap
                break
            end
        end
        gap += diff
    end
    # @info "All equal!"
    return TIequality, gap
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
    subproblems = [Subproblem(t, problem, chosentrucks) for t in 1:batchsize]

    optsolpertruck = [Dict{Symbol, Any}() for t in 1:nbtrucks]
    
    TIequality, gap = are_all_TI_equal(TIvalues, TIbar, nbchosentrucks, eps, retgap=true) # TODO leave possibility to user to enable/disable gap 
    oldgap = gap
    # @debug "TIequality" TIequality
    # alltrucksdone = falses(nbchosentrucks)
    firstpass = true
    iteration = 0
    while (!TIequality && oldgap >= gap) || firstpass 
    # while !TIequality
        firstpass = false
        begin
            # # @info "Convergence gap" sum(abs.(TIvalues[t, :, :] .- TIbar))
            pt = 1
            iteration +=1
            printstyled("Uzawa -- iteration: ", iteration, " -- gap: ", gap, "\n", color=:blue)
            for pt in 1:nbplannedtrucks
                display_progress(pt, nbplannedtrucks, n=20, name="Planned truck")
                # @sync @distributed for b in 0:batchsize-1 # Threads
                for b in 0:batchsize-1 # Threads
                    display_progress(b+1, batchsize, n=20, name="Thread")
                    # @debug "Thread $b" 
                    ### Assign a planned truck (and corresponding extra trucks) to a thread/subproblem instance
                    ptb = pt + b
                    if ptb > nbplannedtrucks
                        continue
                    end
                    chosentrucks_ptb = filter(x -> x in chosentrucks, problem[:truckindices][ptb])
                    for i in chosentrucks_ptb
                        idx = findfirst(==(i), chosentrucks)
                        display_progress(findfirst(==(i), chosentrucks_ptb), length(chosentrucks_ptb), n=20, name="Truck $i")
                        # @debug "truck" i
                        # If the conscensus is that there are no items in this truck, the 
                        # truck agrees and there is no need to solve it.
                        if sum(TIbar[idx, :]) == 0
                        # if sum(TIbar[i, :]) == 0 || sum(problem[:TR][i, :]) == 0
                            # @debug "=> Truck empty"
                            # setvalueTI!(subproblems[t], TIbar)
                            TIvalues[idx, :, :] .= round.(TIbar) # TODO round is risky?
                            # optsolpertruck[i][:S] = nothing
                            # optsolpertruck[i][:SG] = nothing
                            # optsolpertruck[i][:SL] = nothing
                            # optsolpertruck[i][:SW] = nothing
                            # optsolpertruck[i][:SH] = nothing

                        else
                            # @debug "=> Optimizing subproblem..."
                            changeS = !(truck(subproblems[b+1]) in problem[:truckindices][ptb])
                            # @debug begin
                            #     @debug "Changing truck of subproblem $(b+1)"
                            #     @debug "Was previously $(truck(subproblems[b+1]))"
                            #     @debug "Is now $i"
                            #     @debug "ptb" ptb
                            #     @debug "problem[:truckindices][ptb]" problem[:truckindices][ptb]
                            #     @debug "ReBuilding from scratch: $(changeS ? "yes" : "no" )"
                            # end
                            # @debug "Changing truck..."
                            # changetruck!(i, subproblems[b+1], chosentrucks, changeS=changeS)
                            if changeS
                                setmodel!(subproblems[b+1], Model(optimizer(getproblem(subproblems[b+1]))))
                                buildTSImodel!(model(subproblems[b+1]), getproblem(subproblems[b+1]), i, chosentrucks)
                            else
                                changetruck!(i, subproblems[b+1], chosentrucks)
                            end
                            # Solve the deterministic minimization problem for truck t with
                            # a penalization + kappa[t, k] * (TI[t, k+1] - TIbar[k])
                            # and obtain optimal first decision TI[t, k+1]
                            # @debug "Adding penalization..."
                            upd_penalization!(subproblems[b+1], TIbar, kappas[idx])
                            # @debug "Setting start values..."
                            set_start_value.(model(subproblems[b+1])[:TI], TIbar)
                            # @debug "Optimizing..."
                            # @debug begin
                            #     # @debug "TR" problem[:TR]
                            #     # show(model(subproblems[b+1]))
                            #     # print(model(subproblems[b+1]))
                            #     write_to_file(model(subproblems[b+1]), "model.lp")
                            # end
                            # @debug begin
                            # begin
                            #     # error("Debug stop")
                            #     undo = relax_integrality(model(subproblems[b+1]))
                            #     set_optimizer(model(subproblems[b+1]), Clp.Optimizer)
                            #     # cons = find_problematic_constraint!(model(subproblems[b+1]))
                            #     # open("tmp_infeasibilities.txt", "w") do io
                            #     #     show(io, "text/plain", cons)
                            #     # end
                            #     # error("Stop")
                            # end
                            optimize!(model(subproblems[b+1]))
                            # error("Stop") # debug
                            # upd_valueTI!(subproblems[t])
                            TIvalues[idx, :, :] .= round.(value.(model(subproblems[b+1])[:TI])) # TODO round is risky?
                            if haskey(model(subproblems[b+1]), :S)
                                optsolpertruck[i][:S] = copy(value.(model(subproblems[b+1])[:S]))
                            end
                            # optsolpertruck[i][:SG] = copy(value.(model(subproblems[b+1])[:SG]))
                            # optsolpertruck[i][:SL] = copy(value.(model(subproblems[b+1])[:SL]))
                            # optsolpertruck[i][:SW] = copy(value.(model(subproblems[b+1])[:SW]))
                            # optsolpertruck[i][:SH] = copy(value.(model(subproblems[b+1])[:S])) * problem[:IH]
                            if haskey(model(subproblems[b+1]), :SXo)
                                optsolpertruck[i][:SXo] = copy(value.(model(subproblems[b+1])[:SXo]))
                            end
                            if haskey(model(subproblems[b+1]), :SXe)
                                optsolpertruck[i][:SXe] = copy(value.(model(subproblems[b+1])[:SXe]))
                            end
                            if haskey(model(subproblems[b+1]), :SYo)
                                optsolpertruck[i][:SYo] = copy(value.(model(subproblems[b+1])[:SYo]))
                            end
                            if haskey(model(subproblems[b+1]), :SYe)
                                optsolpertruck[i][:SYe] = copy(value.(model(subproblems[b+1])[:SYe]))
                            end
                            if haskey(model(subproblems[b+1]), :SZe)
                                optsolpertruck[i][:SZe] = copy(value.(model(subproblems[b+1])[:SZe]))
                            end
                            if haskey(model(subproblems[b+1]), :SO)
                                optsolpertruck[i][:SO] = copy(value.(model(subproblems[b+1])[:SO]))
                            end
                            if haskey(model(subproblems[b+1]), :SU)
                                optsolpertruck[i][:SU] = copy(value.(model(subproblems[b+1])[:SU]))
                            end
                            if haskey(model(subproblems[b+1]), :SK)
                                optsolpertruck[i][:SK] = copy(value.(model(subproblems[b+1])[:SK]))
                            end
                            optsolpertruck[i][:ST] = i # ? TODO
                        end
                    end
                    clearnlines(1)
                end
                clearnlines(1)
                pt += batchsize
            end
            # clearnlines(1)
            clearnlines(2)
            # printstyled("Done", "\n", color=:green)

            # Update the mean first decisions:
            # TIbar[k+1] = sum([pi[t] * TI[t, k+1] for t in bold_T])
            # TIbar = sum([value.(model(subproblems[t])[:TI]) for t in 1:nbchosentrucks])
            TIbar = sum([TIvalues[t, :, :] for t in 1:nbchosentrucks])/nbchosentrucks

            # Update the multipliers by
            for t in 1:nbchosentrucks
                kappas[t] = kappas[t] .+ delta * sum(TIvalues[t, :, :] .- TIbar)
            end
            oldgap = gap
            TIequality, gap = are_all_TI_equal(TIvalues, TIbar, nbchosentrucks, eps; retgap=true)
        end
    end


    return TIbar, optsolpertruck
end

function changetruck!(t, subproblem::Subproblem, chosentrucks)
    
    return buildTSImodel!(model(subproblem), getproblem(subproblem), t, chosentrucks; replace=true)

    # ## Create model
    # submodel = model(subproblem)
    # nbcandidateitems = max(1, sum(subproblem[:TR][chosentrucks[t], :]))
    # # nbcandidateitems = max([sum(problem[:TR][t2, :]) for t2 in 1:nbtrucks]...)
    # nbstacks = max(1, nbcandidateitems) # in case no candidate item, at least equal to 1
    # nbitems = subproblem[:nbitems]
    # nbtrucks = subproblem[:nbtrucks]
    # nbchosentrucks = length(chosentrucks)
    # nbplannedtrucks = subproblem[:nbplannedtrucks]
    # # nbplants = subproblem[:nbplants]
    # nbplantdocks = subproblem[:nbplantdocks]
    # nbsuppliers = subproblem[:nbsuppliers]
    # nbsupplierdocks = subproblem[:nbsupplierdocks]
    # # nbextratrucks = nbtrucks - nbplannedtrucks
    # nbextratrucks = length(filter(x -> !(x in 1:nbplannedtrucks), chosentrucks))
    # subproblem.t = t
    # # subproblem = Subproblem(t, subproblem, submodel, optimizer, nbstacks, Matrix{Union{Missing, Bool}}(missing, nbtrucks, nbtrucks))

    # if changeS
    #     # # @info "Replacing S..."
    #     delete.(submodel, submodel[:S])
    #     JuMP.unregister.(submodel, :S)
    #     @variable(submodel, S[1:nbstacks, 1:nbcandidateitems], lower_bound = 0, upper_bound = 1, Bin)

    #     # # @info "Replacing SS..."
    #     delete.(submodel, submodel[:SS])
    #     JuMP.unregister.(submodel, :SS)
    #     @variable(submodel, SS[1:nbstacks] >= 0)
    #     # # @info "Replacing SP..."
    #     delete.(submodel, submodel[:SP])
    #     JuMP.unregister.(submodel, :SP)
    #     @variable(submodel, SP[1:nbstacks], lower_bound = 0)
    #     # # @info "Replacing SK..."
    #     delete.(submodel, submodel[:SK])
    #     JuMP.unregister.(submodel, :SK)
    #     @variable(submodel, SK[1:nbstacks, 1:nbsupplierdocks], lower_bound = 0, upper_bound = 1) # No need for Bin because it is constrained to be equal to their items' which are integer
    #     # # @info "Replacing SU..."
    #     delete.(submodel, submodel[:SU])
    #     JuMP.unregister.(submodel, :SU)
    #     @variable(submodel, SU[1:nbstacks, 1:nbsuppliers], lower_bound = 0, upper_bound = 1)
    #     # # @info "Replacing SO..."
    #     delete.(submodel, submodel[:SO])
    #     JuMP.unregister.(submodel, :SO)
    #     @variable(submodel, SO[1:nbstacks] >= 0)
    #     # # @info "Replacing SXe..."
    #     delete.(submodel, submodel[:SXe])
    #     JuMP.unregister.(submodel, :SXe)
    #     @variable(submodel, SXe[1:nbstacks] >= 0)
    #     # # @info "Replacing SXo..."
    #     delete.(submodel, submodel[:SXo])
    #     JuMP.unregister.(submodel, :SXo)
    #     @variable(submodel, SXo[1:nbstacks] >= 0)
    #     # # @info "Replacing SYe..."
    #     delete.(submodel, submodel[:SYe])
    #     JuMP.unregister.(submodel, :SYe)
    #     @variable(submodel, SYe[1:nbstacks] >= 0)
    #     # # @info "Replacing SYo..."
    #     delete.(submodel, submodel[:SYo])
    #     JuMP.unregister.(submodel, :SYo)
    #     @variable(submodel, SYo[1:nbstacks] >= 0)
    #     # # @info "Replacing SZe..."
    #     delete.(submodel, submodel[:SZe])
    #     JuMP.unregister.(submodel, :SZe)
    #     @variable(submodel, SZe[1:nbstacks] >= 0)
    #     # # @info "Replacing betaM..."
    #     delete.(submodel, submodel[:betaM])
    #     JuMP.unregister.(submodel, :betaM)
    #     @variable(submodel, betaM[1:convert(Int, nbstacks*(nbstacks+1)/2)] >= 0)
    #     # # @info "Replacing betaP..."
    #     delete.(submodel, submodel[:betaP])
    #     JuMP.unregister.(submodel, :betaP)
    #     @variable(submodel, betaP[1:convert(Int, nbstacks*(nbstacks+1)/2)] >= 0)
    #     # # @info "Replacing nu..."
    #     delete.(submodel, submodel[:nu])
    #     JuMP.unregister.(submodel, :nu)
    #     @variable(submodel, nu[1:nbstacks] >= 0)
    #     # # @info "Replacing tau..."
    #     delete.(submodel, submodel[:tau])
    #     JuMP.unregister.(submodel, :tau)
    #     @variable(submodel, tau[1:nbstacks] >= 0)
    #     # # @info "Replacing phi..."
    #     delete.(submodel, submodel[:phi])
    #     JuMP.unregister.(submodel, :phi)
    #     @variable(submodel, phi[1:nbstacks] >= 0)
    #     # # @info "Replacing SG..."
    #     delete.(submodel, submodel[:SG])
    #     JuMP.unregister.(submodel, :SG)
    #     @variable(submodel, SG[1:nbstacks, 1:nbplantdocks] >= 0, upper_bound = 1)
    
    #     # # @info "Replacing SL..."
    #     delete.(submodel, submodel[:SL])
    #     JuMP.unregister.(submodel, :SL)
    #     @variable(submodel, SL[1:nbstacks] >= 0)
    #     delete.(submodel, submodel[:SW])
    #     JuMP.unregister.(submodel, :SW)
    #     @variable(submodel, SW[1:nbstacks] >= 0)
    #     # # @info "Replacing Z..."
    #     delete.(submodel, submodel[:Z])
    #     JuMP.unregister.(submodel, :Z)
    #     @variable(submodel, Z[1:nbstacks, 1:nbcandidateitems], lower_bound = 0, upper_bound = 1, Bin)
    #     # # @info "Replacing mu..."
    #     delete.(submodel, submodel[:mu])
    #     JuMP.unregister.(submodel, :mu)
    #     @variable(submodel, mu[1:convert(Int, nbstacks*(nbstacks+1)/2)], lower_bound = 0, upper_bound = 1, Bin)
    #     # # @info "Replacing eta..."
    #     delete.(submodel, submodel[:eta])
    #     JuMP.unregister.(submodel, :eta)
    #     @variable(submodel, eta[1:convert(Int, nbstacks*(nbstacks+1)/2)], lower_bound = 0, upper_bound = 1, Bin)
    #     # # @info "Replacing xi..."
    #     delete.(submodel, submodel[:xi])
    #     JuMP.unregister.(submodel, :xi)
    #     @variable(submodel, xi[1:convert(Int, nbstacks*(nbstacks+1)/2)], lower_bound = 0, upper_bound = 1, Bin)
    #     # # @info "Replacing chi..."
    #     delete.(submodel, submodel[:chi])
    #     JuMP.unregister.(submodel, :chi)
    #     @variable(submodel, chi[1:convert(Int, nbstacks*(nbstacks+1)/2)], lower_bound = 0, upper_bound = 1, Bin)
    #     # # @info "Replacing r..."
    #     delete.(submodel, submodel[:r])
    #     JuMP.unregister.(submodel, :r)
    #     @variable(submodel, r[1:convert(Int, nbstacks*(nbstacks+1)/2)], lower_bound = 0, upper_bound = 1, Bin)
    #     # # @info "Replacing sigma1..."
    #     delete.(submodel, submodel[:sigma1])
    #     JuMP.unregister.(submodel, :sigma1)
    #     @variable(submodel, sigma1[1:convert(Int, nbstacks*(nbstacks+1)/2)], lower_bound = 0, upper_bound = 1, Bin)
    #     # # @info "Replacing sigma2..."
    #     delete.(submodel, submodel[:sigma2])
    #     JuMP.unregister.(submodel, :sigma2)
    #     @variable(submodel, sigma2[1:convert(Int, nbstacks*(nbstacks+1)/2)], lower_bound = 0, upper_bound = 1, Bin)
    #     # # @info "Replacing sigma3..."
    #     delete.(submodel, submodel[:sigma3])
    #     JuMP.unregister.(submodel, :sigma3)
    #     @variable(submodel, sigma3[1:convert(Int, nbstacks*(nbstacks+1)/2)], lower_bound = 0, upper_bound = 1, Bin)
    #     # # @info "Replacing Q..."
    #     delete.(submodel, submodel[:Q])
    #     JuMP.unregister.(submodel, :Q)
    #     @variable(submodel, Q[1:nbstacks, 1:nbsuppliers], lower_bound = 0, Int)
    #     # # @info "Replacing V..."
    #     delete.(submodel, submodel[:V])
    #     JuMP.unregister.(submodel, :V)
    #     @variable(submodel, V[1:nbstacks, 1:nbsupplierdocks], lower_bound = 0, Int)
    #     # # @info "Replacing W..."
    #     delete.(submodel, submodel[:W])
    #     JuMP.unregister.(submodel, :W)
    #     @variable(submodel, W[1:nbstacks, 1:nbplantdocks], lower_bound = 0, Int)
    #     # # @info "Replacing Gl..."
    #     delete.(submodel, submodel[:Gl])
    #     JuMP.unregister.(submodel, :Gl)
    #     @variable(submodel, Gl[1:nbstacks, 1:nbcandidateitems], lower_bound = 0, Int)
    #     # # @info "Replacing Gr..."
    #     delete.(submodel, submodel[:Gr])
    #     JuMP.unregister.(submodel, :Gr)
    #     @variable(submodel, Gr[1:nbstacks, 1:nbcandidateitems], lower_bound = 0, Int)
    #     # # @info "Replacing DL..."    
    #     delete.(submodel, submodel[:DL])
    #     JuMP.unregister.(submodel, :DL)
    #     @variable(submodel, DL[1:nbstacks, 1:nbcandidateitems], lower_bound = 0)
    #     # @debug begin
    #     #     if length(DL) == 0
    #     #         @debug "DL" DL
    #     #         @debug "S" S
    #     #         sleep(20)
    #     #     end
    #     # end
    #     # # @info "Replacing DW..."
    #     delete.(submodel, submodel[:DW])
    #     JuMP.unregister.(submodel, :DW)
    #     @variable(submodel, DW[1:nbstacks, 1:nbcandidateitems], lower_bound = 0)
    
    #     # # @info "Replacing lambda..."
    #     delete.(submodel, submodel[:lambda])
    #     JuMP.unregister.(submodel, :lambda)
    #     @variable(submodel, lambda[1:convert(Int, nbstacks * (nbstacks+1)/2)], lower_bound = 0, Bin)
                   
    # end
    # # # @info "Computing parameters..."
    # MI4 = [[max(subproblem[:TU][:, j]...) for j in 1:first(size(subproblem[:TU][1, :]))] for j in 1:first(size(submodel[:TI][:, 1]))]
    # MI4 = map(Int, (max(subproblem[:TU][:, i]...) for i = 1:first(size(subproblem[:TU][1, :])), j = 1:first(size(submodel[:TI][:, 1]))))

    # MZ = max(subproblem[:IS]...) + 1.0

    # MQ = max(subproblem[:IU]...) + 1.0

    # # MH = max(IP...) + 1.0

    # MV = max(subproblem[:IK]...) + 1.0

    # MW = max(subproblem[:IPD]...) + 1.0

    # # MG = max(skipmissing(subproblem[:_IO])...) + 1.0
    # MG = 5.0

    # MDL =  max(subproblem[:IL]...) + 1.0
    # MDW =  max(subproblem[:IW]...) + 1.0

    # # Meta = nbtrucks * Mtau/10

    # MTL = Matrix{Float64}(undef, nbchosentrucks, 1)
    # MTW = Matrix{Float64}(undef, nbchosentrucks, 1)

    # MTL = max(subproblem[:TL]...) + 1.0
    # MTW = max(subproblem[:TW]...) + 1.0
    # # if true
    # # # if typeof(MTL) != Vector{T} where {T <: Real}
    # #     throw(TypeError(MTL, "MTL must be of type Vector{T} where {T <: Real}", Vector{T} where {T <: Real}, typeof(MTL)))
    # # end

    # MTE = max(skipmissing(subproblem[:TE])...)

    # MTKE = max(subproblem[:TKE]...)

    # MTGE = max(subproblem[:TGE]...)

    # MTW = max(subproblem[:TW]...) + 1.0

    # # MPsi = nbitems

    # # Mlambda = 2.0 * subproblem[:TL][t] + 1.0
    # Mlambda = 2.0 * max(subproblem[:TL]...) + 1.0

    # Mzeta = nbitems

    # # SZo = 0.0

    # # MST = 2.0

    # # MOmega = 2.0

    # Mmu = 2.0

    # # Mtau = 10.0

    # MS = 2.0

    # Xi1 = falses(convert(Int, nbstacks*(nbstacks+1)/2), nbstacks)
    # Xi2 = falses(convert(Int, nbstacks*(nbstacks+1)/2), nbstacks)

    # epsilon = 0.001


    # fillXi1!(Xi1)
    # fillXi2!(Xi2)

    # ## Add constraints
    # # # @info "Replacing constraints..."
    # if haskey(JuMP.object_dictionary(submodel), :cZetaT2)
    #     delete.(submodel, submodel[:cZetaT2])
    #     JuMP.unregister.(submodel, :cZetaT2)
    # end
    # if haskey(JuMP.object_dictionary(submodel), :cZetaT3)
    #     delete.(submodel, submodel[:cZetaT3])
    #     JuMP.unregister.(submodel, :cZetaT3)
    # end
    # if haskey(JuMP.object_dictionary(submodel), :cZetaE2)
    #     delete.(submodel, submodel[:cZetaE2])
    #     JuMP.unregister.(submodel, :cZetaE2)
    # end
    # if haskey(JuMP.object_dictionary(submodel), :cZetaE3)
    #     delete.(submodel, submodel[:cZetaE3])
    #     JuMP.unregister.(submodel, :cZetaE3)
    # end
    # if t <= nbplannedtrucks
    #     # # @info "Replacing cZetaT2..."
    #     @constraint(submodel, cZetaT2, submodel[:zetaT] * Mzeta .>= submodel[:TI][1:nbplannedtrucks, :] * vones(Int8, nbitems))
    #     @constraint(submodel, cZetaT3, -(1 .- submodel[:zetaT]) * Mzeta .+ 1 .<= submodel[:TI][1:nbplannedtrucks, :] * vones(Int8, nbitems)) # TODO correct
    # else
    #     # # @info "Replacing cZetaE2..."
    #     @constraint(submodel, cZetaE2, submodel[:zetaE] * Mzeta .>= submodel[:TI][nbplannedtrucks+1:end, :] * vones(Int8, nbitems))
    #     @constraint(submodel, cZetaE3, -(1 .- submodel[:zetaE]) * Mzeta .+ 1 .<= submodel[:TI][nbplannedtrucks+1:end, :] * vones(Int8, nbitems))
    # end
    # # # @info "Replacing cTI_TR..."
    # delete.(submodel, submodel[:cTI_TR])
    # JuMP.unregister.(submodel, :cTI_TR)
    # @constraint(submodel, cTI_TR, submodel[:TI][t, :] .<= subproblem[:TR][chosentrucks[t], :])

    # icandidates = findall((x) -> x == 1, subproblem[:TR][chosentrucks[t], :])
    # # # @info "Replacing cTI_1_1..."
    # delete.(submodel, submodel[:cTI_1_1])
    # JuMP.unregister.(submodel, :cTI_1_1)
    # @constraint(submodel, cTI_1_1, transpose(submodel[:TI])[filter(x -> x in icandidates, 1:nbitems), :] * vones(Int8, nbchosentrucks) + submodel[:GI][filter(x -> x in icandidates, 1:nbitems)] .<= vones(Int8, size(transpose(submodel[:TI])[filter(x -> x in icandidates, 1:nbitems), :], 1)))
    # # @debug "cTI_1_1" cTI_1_1[icandidates[1], :]
    # # # @info "Replacing cS_TI..."
    # delete.(submodel, submodel[:cS_TI])
    # JuMP.unregister.(submodel, :cS_TI)
    # # @debug begin
    # #     @debug "transpose(submodel[:S])" transpose(submodel[:S])
    # #     @debug "submodel[:TI][chosentrucks[t], filter(x -> x in icandidates, 1:nbitems)]" submodel[:TI][chosentrucks[t], filter(x -> x in icandidates, 1:nbitems)]
    # #     @debug "nbcandidateitems" nbcandidateitems
    # #     @debug "icandidates" icandidates
    # #     @debug "submodel[:S]" submodel[:S]
    # #     @debug "changeS" changeS
    # # end 
    # @constraint(submodel, cS_TI, transpose(submodel[:S]) * vones(Int8, nbstacks) .== submodel[:TI][chosentrucks[t], filter(x -> x in icandidates, 1:nbitems)])

    # if changeS
    #     # # @info "Replacing cZ_S_MZ..."
    #     delete.(submodel, submodel[:cZ_S_MZ])
    #     JuMP.unregister.(submodel, :cZ_S_MZ)
    #     @constraint(submodel, cZ_S_MZ, submodel[:Z] .<= submodel[:S] * MZ)
    # end
    # # # # @info "Replacing cS_TI_MS..."
    # # delete.(submodel, submodel[:cS_TI_MSleft])
    # # JuMP.unregister.(submodel, :cS_TI_MSleft)
    # # @constraint(submodel, cS_TI_MSleft, -(vones(Int8, nbcandidateitems) - submodel[:TI][t, filter(x -> x in icandidates, 1:nbitems)]) * MS .<= transpose(submodel[:S]) * vones(Int8, nbcandidateitems) - vones(Int8, nbcandidateitems))
    # # delete.(submodel, submodel[:cS_TI_MSright])
    # # JuMP.unregister.(submodel, :cS_TI_MSright)
    # # @constraint(submodel, cS_TI_MSright, transpose(submodel[:S]) * vones(Int8, nbcandidateitems) - vones(Int8, nbcandidateitems) .<= (vones(Int8, nbcandidateitems) - submodel[:TI][t, filter(x -> x in icandidates, 1:nbitems)]) * MS)

    # # # @info "Replacing cS_IS_Z..."
    # delete.(submodel, submodel[:cS_IS_Z])
    # JuMP.unregister.(submodel, :cS_IS_Z)
    # # @debug begin
    # #     @debug "submodel[:S]" submodel[:S]
    # #     @debug "subproblem[:IS][filter(x -> x in icandidates, 1:nbitems)]" subproblem[:IS][filter(x -> x in icandidates, 1:nbitems)]
    # # end
    # @constraint(submodel, cS_IS_Z, submodel[:S] * subproblem[:IS][filter(x -> x in icandidates, 1:nbitems)] .== submodel[:Z] * vones(Int8, size(submodel[:Z])[1]))

    # if changeS
    #     # # @info "Replacing cZ_SS_S..."
    #     delete.(submodel, submodel[:cZ_SS_Sleft])
    #     JuMP.unregister.(submodel, :cZ_SS_Sleft)
    #     @constraint(submodel, cZ_SS_Sleft, -MZ * (-submodel[:S] .+ 1) .<= submodel[:Z] .- hcat([submodel[:SS] for i in 1:size(submodel[:Z])[2]]...))
    #     delete.(submodel, submodel[:cZ_SS_Sright])
    #     JuMP.unregister.(submodel, :cZ_SS_Sright)
    #     @constraint(submodel, cZ_SS_Sright, submodel[:Z] .- hcat([submodel[:SS] for i in 1:size(submodel[:Z])[2]]...) .<= MZ * (-submodel[:S] .+ 1))
    #     # # @info "Replacing cQ_SU..."
    #     delete.(submodel, submodel[:cQ_SU])
    #     JuMP.unregister.(submodel, :cQ_SU)
    #     @constraint(submodel, cQ_SU, submodel[:Q] .<= submodel[:SU] * MQ)

    # end

    # # # @info "Replacing cS_IU_Q..."
    # delete.(submodel, submodel[:cS_IU_Q])
    # JuMP.unregister.(submodel, :cS_IU_Q)
    # @constraint(submodel, cS_IU_Q, submodel[:S] * subproblem[:IU][filter(x -> x in icandidates, 1:nbitems)] .== submodel[:Q])
    
    # if changeS
    #     # # @info "Replacing cQ_SU_S..."
    #     delete.(submodel, submodel[:cQ_SU_Sleft])
    #     JuMP.unregister.(submodel, :cQ_SU_Sleft)
    #     @constraint(submodel, cQ_SU_Sleft, -MQ * (1 .- submodel[:SU]) .<= submodel[:Q] .- hcat([submodel[:S] * vones(Int8, size(submodel[:S], 2)) for i in 1:size(submodel[:Q], 2)]...))
    #     delete.(submodel, submodel[:cQ_SU_Sright])
    #     JuMP.unregister.(submodel, :cQ_SU_Sright)
    #     @constraint(submodel, cQ_SU_Sright, submodel[:Q] .- hcat([submodel[:S] * vones(Int8, size(submodel[:S], 2)) for i in 1:size(submodel[:Q], 2)]...) .<= MQ * (1 .- submodel[:SU]))
    #     # # @info "Replacing cV_SK..."
    #     delete.(submodel, submodel[:cV_SK])
    #     JuMP.unregister.(submodel, :cV_SK)
    #     @constraint(submodel, cV_SK, submodel[:V] .<= submodel[:SK] * MV)
        
    #     # # @info "Replacing cS_IK_V..."
    #     delete.(submodel, submodel[:cS_IK_V])
    #     JuMP.unregister.(submodel, :cS_IK_V)
    #     @constraint(submodel, cS_IK_V, submodel[:S] * subproblem[:IK][filter(x -> x in icandidates, 1:nbitems)] .== submodel[:V])
    # end

    # if changeS
    #     # # @info "Replacing cV_SK_S..."
    #     delete.(submodel, submodel[:cV_SK_Sleft])
    #     JuMP.unregister.(submodel, :cV_SK_Sleft)
    #     @constraint(submodel, cV_SK_Sleft, -MV * (-submodel[:SK] .+ 1) .<= submodel[:V] .- hcat([submodel[:S] * vones(Int8, size(submodel[:S], 2)) for i in 1:size(submodel[:V], 2)]...))
    #     delete.(submodel, submodel[:cV_SK_Sright])
    #     JuMP.unregister.(submodel, :cV_SK_Sright)
    #     @constraint(submodel, cV_SK_Sright, submodel[:V] .- hcat([submodel[:S] * vones(Int8, size(submodel[:S], 2)) for i in 1:size(submodel[:V], 2)]...) .<= MV * (-submodel[:SK] .+ 1))
    #     # # @info "Replacing cW_SG..."
    #     delete.(submodel, submodel[:cW_SG])
    #     JuMP.unregister.(submodel, :cW_SG)
    #     @constraint(submodel, cW_SG, submodel[:W] .<= submodel[:SG] * MW)
        
    #     # # @info "Replacing cS_IPD_W..."
    #     delete.(submodel, submodel[:cS_IPD_W])
    #     JuMP.unregister.(submodel, :cS_IPD_W)
    #     @constraint(submodel, cS_IPD_W, submodel[:S] * subproblem[:IPD][filter(x -> x in icandidates, 1:nbitems)] .== submodel[:W])

    #     # # @info "Replacing cW_SG_S..."
    #     delete.(submodel, submodel[:cW_SPD_Sleft])
    #     JuMP.unregister.(submodel, :cW_SPD_Sleft)
    #     @constraint(submodel, cW_SPD_Sleft, -MW * (-submodel[:SG] .+ 1) .<= submodel[:W] .- hcat([submodel[:S] * vones(Int8, size(submodel[:S], 2)) for i in 1:size(submodel[:W], 2)]...))
    #     delete.(submodel, submodel[:cW_SPD_Sright])
    #     JuMP.unregister.(submodel, :cW_SPD_Sright)
    #     @constraint(submodel, cW_SPD_Sright, submodel[:W] .- hcat([submodel[:S] * vones(Int8, size(submodel[:S], 2)) for i in 1:size(submodel[:W], 2)]...) .<= MW * (-submodel[:SG] .+ 1))

    #     # # @info "Replacing cGl_S..."
    #     delete.(submodel, submodel[:cGl_S])
    #     JuMP.unregister.(submodel, :cGl_S)
    #     @constraint(submodel, cGl_S, submodel[:Gl] .<= submodel[:S] * MG)
    #     # # @info "Replacing cGr_S..."
    #     delete.(submodel, submodel[:cGr_S])
    #     JuMP.unregister.(submodel, :cGr_S)
    #     @constraint(submodel, cGr_S, submodel[:Gr] .<= submodel[:S] * MG)

    #     # # @info "Replacing cGl_Gr..."
    #     delete.(submodel, submodel[:cGl_Gr])
    #     JuMP.unregister.(submodel, :cGl_Gr)
    #     @constraint(submodel, cGl_Gr, submodel[:Gl] .== submodel[:Gr])

    #     # # @info "Replacing cGr_SO_S..."
    #     delete.(submodel, submodel[:cGr_SO_Sleft])
    #     JuMP.unregister.(submodel, :cGr_SO_Sleft)
    #     @constraint(submodel, cGr_SO_Sleft, -MG * (-submodel[:S] .+ 1) .<= submodel[:Gr] .- hcat([submodel[:SO] for i in 1:nbcandidateitems]...))
    #     delete.(submodel, submodel[:cGr_SO_Sright])
    #     JuMP.unregister.(submodel, :cGr_SO_Sright)
    #     @constraint(submodel, cGr_SO_Sright, submodel[:Gr] .- hcat([submodel[:SO] for i in 1:nbcandidateitems]...) .<= MG * (-submodel[:S] .+ 1))
        
    #     # # @info "Replacing cGl_IOV_S..."
    #     delete.(submodel, submodel[:cGl_IOV_Sleft])
    #     JuMP.unregister.(submodel, :cGl_IOV_Sleft)
    #     @constraint(submodel, cGl_IOV_Sleft, -MG * (-submodel[:S] .+ 1) .<= submodel[:Gl] .- vcat([transpose(submodel[:IOV][filter(x -> x in icandidates, 1:nbitems)]) for i in 1:nbstacks]...))
    #     delete.(submodel, submodel[:cGl_IOV_Sright])
    #     JuMP.unregister.(submodel, :cGl_IOV_Sright)
    #     @constraint(submodel, cGl_IOV_Sright, submodel[:Gl] .- vcat([transpose(submodel[:IOV][filter(x -> x in icandidates, 1:nbitems)]) for i in 1:nbstacks]...) .<= MG * (-submodel[:S] .+ 1))

    #     # # @info "Replacing cDL_S..."
    #     delete.(submodel, submodel[:cDL_S])
    #     JuMP.unregister.(submodel, :cDL_S)
    #     @constraint(submodel, cDL_S, submodel[:DL] .<= submodel[:S] * MDL)
    #     # # @info "Replacing cDL_SL..."
    #     delete.(submodel, submodel[:cDL_SLleft])
    #     JuMP.unregister.(submodel, :cDL_SLleft)
    #     @constraint(submodel, cDL_SLleft, -MDL * (-submodel[:S] .+ 1) .<= submodel[:DL] - hcat([submodel[:SL] for i in 1:size(submodel[:DL])[2]]...))
    #     delete.(submodel, submodel[:cDL_SLright])
    #     JuMP.unregister.(submodel, :cDL_SLright)
    #     # @debug begin
    #     #     @debug "submodel[:DL]" submodel[:DL]
    #     #     @debug "hcat([submodel[:SL] for i in 1:size(submodel[:DL])[2]]...)" hcat([submodel[:SL] for i in 1:size(submodel[:DL])[2]]...)
    #     #     @debug "MDL" MDL
    #     #     @debug "(-submodel[:S] .+ 1)" (-submodel[:S] .+ 1)
    #     # end
    #     @constraint(submodel, cDL_SLright, submodel[:DL] - hcat([submodel[:SL] for i in 1:size(submodel[:DL])[2]]...) .<= MDL * (-submodel[:S] .+ 1))
    #     # # @info "Replacing cDL_S_IL..."
    #     delete.(submodel, submodel[:cDL_S_IL])
    #     JuMP.unregister.(submodel, :cDL_S_IL)
    #     @constraint(submodel, cDL_S_IL, submodel[:DL] * vones(Int8, size(submodel[:S], 1)) .== submodel[:S] * subproblem[:IL][filter(x -> x in icandidates, 1:nbitems)])

    #     # # @info "Replacing cDW_S..."
    #     delete.(submodel, submodel[:cDW_S])
    #     JuMP.unregister.(submodel, :cDW_S)
    #     @constraint(submodel, cDW_S, submodel[:DW] .<= submodel[:S] * MDW)
    #     # # @info "Replacing cDW_SW..."
    #     delete.(submodel, submodel[:cDW_SWleft])
    #     JuMP.unregister.(submodel, :cDW_SWleft)
    #     @constraint(submodel, cDW_SWleft, -MDW * (-submodel[:S] .+ 1) .<= submodel[:DW] - hcat([submodel[:SW] for i in 1:size(submodel[:DW])[2]]...))
    #     delete.(submodel, submodel[:cDW_SWright])
    #     JuMP.unregister.(submodel, :cDW_SWright)
    #     @constraint(submodel, cDW_SWright, submodel[:DW] - hcat([submodel[:SW] for i in 1:size(submodel[:DW])[2]]...) .<= MDW * (-submodel[:S] .+ 1))
        
        
    #     # # @info "Replacing cDW_S_IW..."
    #     delete.(submodel, submodel[:cDW_S_IW])
    #     JuMP.unregister.(submodel, :cDW_S_IW)
    #     @constraint(submodel, cDW_S_IW, submodel[:DW] * vones(Int8, size(submodel[:S], 1)) .== submodel[:S] * subproblem[:IW][filter(x -> x in icandidates, 1:nbitems)])
        
    #     # # @info "Replacing cSZe_S_IH..."
    #     delete.(submodel, submodel[:cSZe_S_IH])
    #     JuMP.unregister.(submodel, :cSZe_S_IH)
    #     @constraint(submodel, cSZe_S_IH, submodel[:SZe] .== submodel[:S]* subproblem[:IH][filter(x -> x in icandidates, 1:nbitems)])

    #     # # @info "Replacing cSXe_SXo_SL_SO..."
    #     delete.(submodel, submodel[:cSXe_SXo_SL_SOleft])
    #     JuMP.unregister.(submodel, :cSXe_SXo_SL_SOleft)
    #     @constraint(submodel, cSXe_SXo_SL_SOleft, submodel[:SXe] - submodel[:SXo] - submodel[:SL] .<= submodel[:SO] * MTL)
    #     delete.(submodel, submodel[:cSXe_SXo_SL_SOright])
    #     JuMP.unregister.(submodel, :cSXe_SXo_SL_SOright)
    #     @constraint(submodel, cSXe_SXo_SL_SOright, submodel[:SXe] - submodel[:SXo] - submodel[:SL] .>= -submodel[:SO] * MTL)
    #     # # @info "Replacing cSYe_SYo_SW_SO..."
    #     delete.(submodel, submodel[:cSYe_SYo_SW_SOleft])
    #     JuMP.unregister.(submodel, :cSYe_SYo_SW_SOleft)
    #     @constraint(submodel, cSYe_SYo_SW_SOleft, submodel[:SYe] - submodel[:SYo] - submodel[:SW] .<= submodel[:SO] * MTW)
    #     delete.(submodel, submodel[:cSYe_SYo_SW_SOright])
    #     JuMP.unregister.(submodel, :cSYe_SYo_SW_SOright)
    #     @constraint(submodel, cSYe_SYo_SW_SOright, submodel[:SYe] - submodel[:SYo] - submodel[:SW] .>= -submodel[:SO] * MTW)
    #     # # @info "Replacing cSXe_SXo_SW_SO..."
    #     delete.(submodel, submodel[:cSXe_SXo_SW_SOleft])
    #     JuMP.unregister.(submodel, :cSXe_SXo_SW_SOleft)
    #     @constraint(submodel, cSXe_SXo_SW_SOleft, submodel[:SXe] - submodel[:SXo] - submodel[:SW] .<= (-submodel[:SO] .+ 1) * MTW)
    #     delete.(submodel, submodel[:cSXe_SXo_SW_SOright])
    #     JuMP.unregister.(submodel, :cSXe_SXo_SW_SOright)
    #     @constraint(submodel, cSXe_SXo_SW_SOright, submodel[:SXe] - submodel[:SXo] - submodel[:SW] .>= -(-submodel[:SO] .+ 1) * MTW)
    #     # # @info "Replacing cSYe_SYo_SL_SO..."
    #     delete.(submodel, submodel[:cSYe_SYo_SL_SOleft])
    #     JuMP.unregister.(submodel, :cSYe_SYo_SL_SOleft)
    #     @constraint(submodel, cSYe_SYo_SL_SOleft, submodel[:SYe] - submodel[:SYo] - submodel[:SL] .<= (-submodel[:SO] .+ 1) * MTL)
    #     delete.(submodel, submodel[:cSYe_SYo_SL_SOright])
    #     JuMP.unregister.(submodel, :cSYe_SYo_SL_SOright)
    #     @constraint(submodel, cSYe_SYo_SL_SOright, submodel[:SYe] - submodel[:SYo] - submodel[:SL] .>= -(-submodel[:SO] .+ 1) * MTL)
    # end

    # # # @info "Replacing cSXe_ST_TL..."
    # delete.(submodel, submodel[:cSXe_ST_TL])
    # JuMP.unregister.(submodel, :cSXe_ST_TL)
    # @constraint(submodel, cSXe_ST_TL, submodel[:SXe] .<= subproblem[:TL][chosentrucks[t]])
    # # # @info "Replacing cSYe_ST_TW..."
    # delete.(submodel, submodel[:cSYe_ST_TW])
    # JuMP.unregister.(submodel, :cSYe_ST_TW)
    # @constraint(submodel, cSYe_ST_TW, submodel[:SYe] .<= subproblem[:TW][chosentrucks[t]])
    # # # @info "Replacing cSZe_ST_TH..."
    # delete.(submodel, submodel[:cSZe_ST_TH])
    # JuMP.unregister.(submodel, :cSZe_ST_TH)
    # @constraint(submodel, cSZe_ST_TH, submodel[:SZe] .<= subproblem[:TH][chosentrucks[t]])
    # if changeS
    #     # # @info "Replacing cSXo_SXo..."
    #     delete.(submodel, submodel[:cSXo_SXo])
    #     JuMP.unregister.(submodel, :cSXo_SXo)
    #     @constraint(submodel, cSXo_SXo, submodel[:SXo][1:end-1] .<= submodel[:SXo][2:end])

    #     # # @info "Adding cXi2SXo_Xi1SXe_betaM_betaP..."
    #     # # @info "Replacing cSXo_SXo..."
    #     delete.(submodel, submodel[:cXi2SXo_Xi1SXe_betaM_betaP])
    #     JuMP.unregister.(submodel, :cXi2SXo_Xi1SXe_betaM_betaP)
    #     @constraint(submodel, cXi2SXo_Xi1SXe_betaM_betaP, (Xi2 * submodel[:SXo]) - (Xi1 * submodel[:SXe]) - submodel[:betaM] + submodel[:betaP] .== -epsilon * vones(Float64, size(Xi1, 1)))
        
        
        
    #     # # @info "Replacing cbetaM_lambda..."
    #     delete.(submodel, submodel[:cbetaM_lambda])
    #     JuMP.unregister.(submodel, :cbetaM_lambda)
    #     @constraint(submodel, cbetaM_lambda, submodel[:betaM] .<= submodel[:lambda] .* Mlambda)
        
    #     # # @info "Replacing cbetaP_lambda..."
    #     delete.(submodel, submodel[:cbetaP_lambda])
    #     JuMP.unregister.(submodel, :cbetaP_lambda)
    #     @constraint(submodel, cbetaP_lambda, submodel[:betaP] .<= (-submodel[:lambda] .+ 1)*Mlambda)
        
    #     # # @info "Replacing cmu_betaM..."
    #     delete.(submodel, submodel[:cmu_betaM])
    #     JuMP.unregister.(submodel, :cmu_betaM)
    #     @constraint(submodel, cmu_betaM, (-submodel[:mu] .+ 1) .<= submodel[:betaM] * Mmu)
        
    #     # # @info "Replacing cXi1SYe_Xi2SYo..."
    #     delete.(submodel, submodel[:cXi1SYe_Xi2SYo])
    #     JuMP.unregister.(submodel, :cXi1SYe_Xi2SYo)
    #     @constraint(submodel, cXi1SYe_Xi2SYo, Xi1 * submodel[:SYe] .<= Xi2 * submodel[:SYo] + submodel[:xi] * MTW + (-submodel[:mu] .+ 1) * MTW)
        
    #     # # @info "Replacing cXi2SYe_Xi1SYo..."
    #     delete.(submodel, submodel[:cXi2SYe_Xi1SYo])
    #     JuMP.unregister.(submodel, :cXi2SYe_Xi1SYo)
    #     @constraint(submodel, cXi2SYe_Xi1SYo, Xi2 * submodel[:SYe] .<= Xi1 * submodel[:SYo] + (-submodel[:xi] .+ 1) * MTW + (-submodel[:mu] .+ 1) * MTW)
        
    # end
    
    
    
    
    # # # @info "Replacing cXi1SU_Xi2SU..."
    # notmissingTE = filter(x -> !ismissing(subproblem[:TE][chosentrucks[t], x]), 1:nbsuppliers)
    # delete.(submodel, submodel[:cXi1SU_Xi2SU])
    # JuMP.unregister.(submodel, :cXi1SU_Xi2SU)
    # @constraint(submodel, cXi1SU_Xi2SU, Xi1 * submodel[:SU][:, notmissingTE] * subproblem[:TE][chosentrucks[t], notmissingTE] .<= Xi2 * submodel[:SU][:, notmissingTE] * subproblem[:TE][chosentrucks[t], notmissingTE])

    # # # @info "Replacing cXi1SU_Xi2SU_chi..."
    # delete.(submodel, submodel[:cXi1SU_Xi2SU_chi])
    # JuMP.unregister.(submodel, :cXi1SU_Xi2SU_chi)
    # @constraint(submodel, cXi1SU_Xi2SU_chi, Xi1*submodel[:SU] - Xi2*submodel[:SU] .>= submodel[:chi] * epsilon - submodel[:r]*MTE - (-submodel[:sigma1] .+ 1) * MTE)

    # # # @info "Replacing cXi2SU_Xi1SU_chi..."
    # delete.(submodel, submodel[:cXi2SU_Xi1SU_chi])
    # JuMP.unregister.(submodel, :cXi2SU_Xi1SU_chi)
    # @constraint(submodel, cXi2SU_Xi1SU_chi, Xi2*submodel[:SU] - Xi1*submodel[:SU] .>= (-submodel[:chi] .+ 1) * epsilon - submodel[:r]*MTE - (-submodel[:sigma1] .+ 1) * MTE)


    # # # @info "Replacing cXi2SK_Xi1SK..."
    # delete.(submodel, submodel[:cXi2SK_Xi1SK])
    # JuMP.unregister.(submodel, :cXi2SK_Xi1SK)
    # notmissingTKE = filter(x -> !ismissing(subproblem[:TKE][chosentrucks[t], x]), 1:nbsupplierdocks)
    # @constraint(submodel, cXi2SK_Xi1SK, Xi2*submodel[:SK][:, notmissingTKE]*subproblem[:TKE][chosentrucks[t], notmissingTKE] .>= Xi1*submodel[:SK][:, notmissingTKE]*subproblem[:TKE][chosentrucks[t], notmissingTKE] - (-submodel[:r] .+ 1) * MTKE)

    # # # @info "Replacing cXi1SK_Xi2SK_chi..."
    # delete.(submodel, submodel[:cXi1SK_Xi2SK_chi])
    # JuMP.unregister.(submodel, :cXi1SK_Xi2SK_chi)
    # @constraint(submodel, cXi1SK_Xi2SK_chi, Xi1*submodel[:SK][:, notmissingTKE]*subproblem[:TKE][chosentrucks[t], notmissingTKE] - Xi2*submodel[:SK][:, notmissingTKE]*subproblem[:TKE][chosentrucks[t], notmissingTKE] .>= submodel[:chi]*epsilon - (-submodel[:sigma2] .+ 1)*MTKE)
    # # # @info "Replacing cXi2SK_Xi1SK_chi..."
    # delete.(submodel, submodel[:cXi2SK_Xi1SK_chi])
    # JuMP.unregister.(submodel, :cXi2SK_Xi1SK_chi)
    # @constraint(submodel, cXi2SK_Xi1SK_chi, Xi2*submodel[:SK][:, notmissingTKE]*subproblem[:TKE][chosentrucks[t], notmissingTKE] - Xi1*submodel[:SK][:, notmissingTKE]*subproblem[:TKE][chosentrucks[t], notmissingTKE] .>= (-submodel[:chi] .+ 1)*epsilon - (-submodel[:sigma2] .+ 1)*MTKE)

    # # # @info "Replacing cXi2SG_Xi1SG..."
    # delete.(submodel, submodel[:cXi2SG_Xi1SG])
    # JuMP.unregister.(submodel, :cXi2SG_Xi1SG)
    # notmissingTGE = filter(x -> !ismissing(subproblem[:TGE][chosentrucks[t], x]), 1:nbplantdocks)
    # @constraint(submodel, cXi2SG_Xi1SG, Xi2*submodel[:SG][:, notmissingTGE]*subproblem[:TGE][chosentrucks[t], notmissingTGE] .>= Xi1*submodel[:SG][:, notmissingTGE]*subproblem[:TGE][chosentrucks[t], notmissingTGE] - (-submodel[:sigma3] .+ 1) * MTGE)

    # # # @info "Replacing csigma1_sigma2_sigma3..."
    # delete.(submodel, submodel[:csigma1_sigma2_sigma3])
    # JuMP.unregister.(submodel, :csigma1_sigma2_sigma3)
    # @constraint(submodel, csigma1_sigma2_sigma3, submodel[:sigma1] + submodel[:sigma2] + submodel[:sigma3] .>= 1)


end
