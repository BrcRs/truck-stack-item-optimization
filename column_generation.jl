using LinearAlgebra

function columngeneration(solvefun!, problem, args...)

    nbplannedtrucks = problem[:nbplannedtrucks]
    truckindices = problem[:truckindices]
    nbitems = problem[:nbitems]
    nbtrucks = problem[:nbtrucks]
    deadtrucks = []
    @debug "truckindices" truckindices
    # @debug begin
    #     sleep(10)
    # end
    # @debug begin
    #     # [truckindices[t][j] for t in 1:nbplannedtrucks for j in filter(x -> truckindices[t][x] in 1:length(truckindices[t]), 1:2)]
    #     for t in 1:nbplannedtrucks
    #         for j in filter(x -> x in 1:length(truckindices[t]), 1:2)
    #             println(t, " ", j)
    #             println(truckindices[t][j])
    #             println()
    #         end
    #     end

    # end
    chosentrucks = [truckindices[t][j] for t in 1:nbplannedtrucks for j in filter(x -> x in 1:length(truckindices[t]), 1:2)]
    firstextratrucks = [truckindices[t][end] + 1 for t in 1:nbplannedtrucks]
    TIbar, optsolpertruck = solvefun!(problem, args..., chosentrucks)

    optsol = Dict{Any, Any}(:TI => copy(TIbar))
    # @debug begin
    #     @debug "[optsol[:TI][t, :] for t in nbplannedtrucks+1:length(chosentrucks)]" [optsol[:TI][t, :] for t in nbplannedtrucks+1:length(chosentrucks)]
    #     @debug "[sum(optsol[:TI][t, :]) for t in nbplannedtrucks+1:length(chosentrucks)]" [sum(optsol[:TI][t, :]) for t in nbplannedtrucks+1:length(chosentrucks)]
    #     @debug "[min(sum(optsol[:TI][t, :]), 1) for t in nbplannedtrucks+1:length(chosentrucks)]" [min(sum(optsol[:TI][t, :]), 1) for t in nbplannedtrucks+1:length(chosentrucks)]
    #     @debug "sum([min(sum(optsol[:TI][t, :]), 1) for t in nbplannedtrucks+1:length(chosentrucks)])" sum([min(sum(optsol[:TI][t, :]), 1) for t in nbplannedtrucks+1:length(chosentrucks)])
    #     @debug "problem[:costextratruck] * sum([min(sum(optsol[:TI][t, :]), 1) for t in nbplannedtrucks+1:length(chosentrucks)])" problem[:costextratruck] * sum([min(sum(optsol[:TI][t, :]), 1) for t in nbplannedtrucks+1:length(chosentrucks)])
    #     @debug "sum(problem[:costextratruck] * [min(sum(optsol[:TI][t, :]), 1) for t in nbplannedtrucks+1:length(chosentrucks)])" sum(problem[:costextratruck] * [min(sum(optsol[:TI][t, :]), 1) for t in nbplannedtrucks+1:length(chosentrucks)])
    # end
    optsol[:value] = problem[:costtransportation] * reduce(+, [min(sum(optsol[:TI][t, :]), 1) for t in 1:nbplannedtrucks], init=0.0) + 
    problem[:costextratruck] * reduce(+, [min(sum(optsol[:TI][t, :]), 1) for t in nbplannedtrucks+1:length(chosentrucks)], init=0.0) + 
    problem[:costinventory] * sum(problem[:IDL] - transpose(optsol[:TI]) * problem[:TDA][filter(x -> x in chosentrucks, 1:nbtrucks)])

    # optsol[:value] = sum(problem[:costtransportation] * [min(sum(optsol[:TI][t, :]), 1) for t in 1:nbplannedtrucks])
    # optsol[:value] = sum(problem[:costextratruck] * [min(sum(optsol[:TI][t, :]), 1) for t in nbplannedtrucks+1:length(chosentrucks)])
    # optsol[:value] = problem[:costinventory] * sum(problem[:IDL] - transpose(optsol[:TI]) * problem[:TDA])
    

    improvement = true

    oldTI = copy(optsol[:TI])
    oldvalue = optsol[:value] * 2


    # while transpose(optsol[:TI]) * vones(Int8, nbtrucks) != vones(Int8, nbitems) || improvement
    while sum(optsol[:TI]) != nbitems || improvement
        printstyled("Column Generation\n", color=:light_green)
        printstyled("Value change: \n", color=:light_green)
        printstyled(optsol[:value] - oldvalue, "\n")
        printstyled("Number of unaffected items: \n", color=:light_green)
        printstyled(nbitems - sum(optsol[:TI]), "\n")
        printstyled("Solving with ", length(chosentrucks), " trucks\n", color=:light_blue)
        # @debug "sort(chosentrucks)" sort(chosentrucks)
        # @debug "optsol[:TI]" optsol[:TI]
        # Upd chosen trucks
        # Add extra trucks for each filled planned trucks
        for (i, t) in enumerate(chosentrucks)
            if sum(optsol[:TI][i, :]) > 0
                if t in 1:nbplannedtrucks
                    if 2 in keys(truckindices[t]) && !(truckindices[t][2] in chosentrucks) && !(truckindices[t][2] in deadtrucks)
                        push!(chosentrucks, truckindices[t][2])
                    end
                    continue
                else
                    if !(t+1 in chosentrucks) && !(t+1 in firstextratrucks) && !(t+1 in deadtrucks)
                        push!(chosentrucks, t+1)
                    end
                    continue
                end
            end
        end


        TIbar, optsolpertruck = solvefun!(problem, args..., chosentrucks)
        optsol[:TI] = copy(TIbar)
        optsol[:variables] = copy(optsolpertruck)

        # @debug begin
        #     @debug "problem[:TDA]" problem[:TDA]
        #     @debug "transpose(optsol[:TI])" transpose(optsol[:TI])
        # end

        optsol[:value] = problem[:costtransportation] * reduce(+, [min(sum(optsol[:TI][t, :]), 1) for t in 1:nbplannedtrucks], init=0.0) + 
        problem[:costextratruck] * reduce(+, [min(sum(optsol[:TI][t, :]), 1) for t in nbplannedtrucks+1:length(chosentrucks)], init=0.0) + 
        problem[:costinventory] * sum(problem[:IDL] - transpose(optsol[:TI]) * problem[:TDA][chosentrucks])
        
       
        # If no unaffected item:
        if sum(optsol[:TI]) == nbitems
            # remove empty unused useless trucks
            for (i, t) in enumerate(chosentrucks)
                if sum(optsol[:TI][i, :]) == 0 && !(t in 1:nbplannedtrucks)
                    push!(deadtrucks, t)
                end
            end
            filter!(x -> !(x in deadtrucks), chosentrucks)
        end

        improvement = optsol[:value] < oldvalue || optsol[:TI] != oldTI
        @debug begin
            @debug "improvement" improvement
            @debug "optsol[:value]" optsol[:value]
            @debug "oldvalue" oldvalue
            @debug "optsol[:TI] != oldTI" optsol[:TI] != oldTI
            sleep(10)
        end
        oldTI = copy(optsol[:TI])
        oldvalue = optsol[:value]


    end

    return optsol

end