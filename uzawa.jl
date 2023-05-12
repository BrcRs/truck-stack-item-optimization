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

"""
    upd_penalization!(subpb::Subproblem, TIbar, kappa)

Update the uzawa penalization `kappa` in the objective function for `subpb`. 
"""
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
    # a big M constraint associated to GI
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

    # Cost of used planned trucks
    @expression(submodel, obj1, sum(subpb[:costtransportation] * submodel[:zetaT]))
    # Cost of used extra trucks
    @expression(submodel, obj2, sum(subpb[:costextratruck] * submodel[:zetaE]))
    # inventory cost
    @expression(submodel, obj3, subpb[:costinventory] * sum(subpb[:IDL] - submodel[:TI][t, :] * subpb[:TDA][t]))
    # uzawa penalization
    @expression(submodel, obj4, kappa * sum(submodel[:TI] - TIbar))
    # penalization of GI
    @expression(submodel, obj5, MGI * sum(submodel[:GI]))

    # @objective(submodel, Min, 
    # sum(subpb[:costtransportation] * submodel[:zetaT]) + 
    # sum(subpb[:costextratruck] * submodel[:zetaE]) + 
    # subpb[:costinventory] * sum(subpb[:IDL] - submodel[:TI][t, :] * subpb[:TDA][t]) + 
    # kappa * sum(submodel[:TI] - TIbar))

    @objective(submodel, Min, obj1 + obj2 + obj3 + obj4 + obj5)

end

"""
    are_all_TI_equal(TIvalues, TIbar, nbtrucks, eps; retgap=false)

Return `true` if the sum of distances between actual TIs of each subproblem and `TIbar`
exceeds `eps`.

If retgap=true, also return the sum of the differences. This considerably slows the
function down since when we don't need the gap value we can return early.
"""
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

"""
    solve_uzawa!(problem::TSIProblem, step::Real, eps::Integer, batchsize::Integer, chosentrucks)

Solve the truck stack item affectation `problem` with an algorithm "à la Uzawa", with trucks chosen in
`chosentrucks`. To each truck is assigned a subproblem restricted to constraints 
and items relevant only to the truck. The variables linking all subproblems are the truck item matrix `TI`. 
All `TI` are considered equal when the sum of the distances doesn't exceed `eps`.
`batchsize` is the number of threads used to  simultaneously solve subproblems, not 
implemented for now TODO.
"""
function solve_uzawa!(problem::TSIProblem, step::Real, eps::Integer, batchsize::Integer, chosentrucks)
    # number of trucks when we count all planned trucks + all possible extra trucks
    nbtrucks = problem[:nbtrucks]
    nbplannedtrucks = problem[:nbplannedtrucks]
    nbitems = problem[:nbitems]
    # nbchosentrucks << nbtrucks
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
            # The idea here is to solve for each planned trucks all similar subproblems, 
            # i.e. a planned trucks and related extra trucks by the same thread as to minimize 
            # the number of models to rebuild 
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
                        # index from global problem (with all extra trucks) will differ from
                        # the ones of the problem based on chosentrucks
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
                            # If the new truck is not of the same kind i.e. 
                            # doesn't have the same data as the previous truck 
                            # of the subproblem the model needs to be rebuilt
                            # (and stack (S) related constraints need to be replaced)
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
                                # Only change truck related constraints and not 
                                # stack constraints if same kind of truck
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
                kappas[t] = kappas[t] .+ step * sum(TIvalues[t, :, :] .- TIbar)
            end
            oldgap = gap
            TIequality, gap = are_all_TI_equal(TIvalues, TIbar, nbchosentrucks, eps; retgap=true)
        end
    end


    return TIbar, optsolpertruck
end

function changetruck!(t, subproblem::Subproblem, chosentrucks)
    
    return buildTSImodel!(model(subproblem), getproblem(subproblem), t, chosentrucks; replace=true)

end
