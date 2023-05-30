using LinearAlgebra

"""
    columngeneration(solvefun!, problem, args...; eps=0.1)

Solve the truck-stack-item problem `problem` thanks to function `solvefun!` by 
iteratively adding trucks to the variables in order to find the best solution 
while limiting the size of the problem.
Return solution if current value is no better than the previous value + eps.
"""
function columngeneration(solvefun!, problem, args...; eps=0.1)

    nbplannedtrucks = problem[:nbplannedtrucks]
    truckindices = problem[:truckindices]
    nbitems = problem[:nbitems]
    nbtrucks = problem[:nbtrucks]

    # trucks to remove from problem
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

    # trucks chosen in the problem
    chosentrucks = [truckindices[t][j] for t in 1:nbplannedtrucks for j in filter(x -> x in 1:length(truckindices[t]), 1:2)]
    firstextratrucks = [truckindices[t][end] + 1 for t in 1:nbplannedtrucks]

    printstyled("Column Generation\n", color=:light_green)
    printstyled("Solving with ", length(chosentrucks), " trucks\n", color=:light_blue)

    # retrieve intial solution via solvefun!, a function of your choice
    TIbar, optsolpertruck = solvefun!(problem, args..., chosentrucks)

    # optsol will store the optimal truck-item (TI) assignment for the global
    # problem, as well as particular solutions for each subproblem (one for each
    # truck)
    optsol = Dict{Any, Any}(:TI => copy(TIbar))
    # @debug begin
    #     @debug "[optsol[:TI][t, :] for t in nbplannedtrucks+1:length(chosentrucks)]" [optsol[:TI][t, :] for t in nbplannedtrucks+1:length(chosentrucks)]
    #     @debug "[sum(optsol[:TI][t, :]) for t in nbplannedtrucks+1:length(chosentrucks)]" [sum(optsol[:TI][t, :]) for t in nbplannedtrucks+1:length(chosentrucks)]
    #     @debug "[min(sum(optsol[:TI][t, :]), 1) for t in nbplannedtrucks+1:length(chosentrucks)]" [min(sum(optsol[:TI][t, :]), 1) for t in nbplannedtrucks+1:length(chosentrucks)]
    #     @debug "sum([min(sum(optsol[:TI][t, :]), 1) for t in nbplannedtrucks+1:length(chosentrucks)])" sum([min(sum(optsol[:TI][t, :]), 1) for t in nbplannedtrucks+1:length(chosentrucks)])
    #     @debug "problem[:costextratruck] * sum([min(sum(optsol[:TI][t, :]), 1) for t in nbplannedtrucks+1:length(chosentrucks)])" problem[:costextratruck] * sum([min(sum(optsol[:TI][t, :]), 1) for t in nbplannedtrucks+1:length(chosentrucks)])
    #     @debug "sum(problem[:costextratruck] * [min(sum(optsol[:TI][t, :]), 1) for t in nbplannedtrucks+1:length(chosentrucks)])" sum(problem[:costextratruck] * [min(sum(optsol[:TI][t, :]), 1) for t in nbplannedtrucks+1:length(chosentrucks)])
    # end
    # Compute value of the current solution
    # This line assumes every planned trucks are always included, but removing 
    # some of them could lead to better performance TODO
    optsol[:value] = problem[:costtransportation] * reduce(+, [min(sum(optsol[:TI][t, :]), 1) for t in 1:nbplannedtrucks], init=0.0) + 
    problem[:costextratruck] * reduce(+, [min(sum(optsol[:TI][t, :]), 1) for t in nbplannedtrucks+1:length(chosentrucks)], init=0.0) + 
    problem[:costinventory] * sum(problem[:IDL] - transpose(optsol[:TI]) * problem[:TDA][filter(x -> x in chosentrucks, 1:nbtrucks)])

    # optsol[:value] = sum(problem[:costtransportation] * [min(sum(optsol[:TI][t, :]), 1) for t in 1:nbplannedtrucks])
    # optsol[:value] = sum(problem[:costextratruck] * [min(sum(optsol[:TI][t, :]), 1) for t in nbplannedtrucks+1:length(chosentrucks)])
    # optsol[:value] = problem[:costinventory] * sum(problem[:IDL] - transpose(optsol[:TI]) * problem[:TDA])
    
    # stores wether the solution of the current iteration is better than the last one
    improvement = true

    oldTI = copy(optsol[:TI])
    oldvalue = optsol[:value] * 2


    # while transpose(optsol[:TI]) * vones(Int8, nbtrucks) != vones(Int8, nbitems) || improvement

    # Iterate until all objects are assigned and there is no more improvement made
    while sum(optsol[:TI]) != nbitems || improvement
        printstyled("Column Generation\n", color=:light_green)
        printstyled("Value change: \n", color=:light_green)
        printstyled(optsol[:value] - oldvalue, "\n")
        printstyled("Number of unassigned items: \n", color=:light_green)
        printstyled(nbitems - sum(optsol[:TI]), "\n")
        printstyled("Solving with ", length(chosentrucks), " trucks ($nbitems items)\n", color=:light_blue)
        oldTI = copy(optsol[:TI])
        oldvalue = optsol[:value]
        # @debug "sort(chosentrucks)" sort(chosentrucks)
        # @debug "optsol[:TI]" optsol[:TI]
        
        # Update chosen trucks
        # Add extra trucks for each filled planned truck
        for (i, t) in enumerate(chosentrucks)
            # if the truck t is non-empty
            if sum(optsol[:TI][i, :]) > 0
                if t in 1:nbplannedtrucks
                    # Add an extra truck if it exists and is not eliminated in deadtrucks
                    if 2 in keys(truckindices[t]) && !(truckindices[t][2] in chosentrucks) && !(truckindices[t][2] in deadtrucks)
                        push!(chosentrucks, truckindices[t][2])
                    end
                else
                    # Same but t is an extra truck
                    if !(t+1 in chosentrucks) && !(t+1 in firstextratrucks) && !(t+1 in deadtrucks)
                        push!(chosentrucks, t+1)
                    end
                end
            end
        end

        # solve subproblems with chosen trucks
        TIbar, optsolpertruck = solvefun!(problem, args..., chosentrucks)
        optsol[:TI] = round.(copy(TIbar))
        optsol[:variables] = copy(optsolpertruck)

        # @debug begin
        #     @debug "problem[:TDA]" problem[:TDA]
        #     @debug "transpose(optsol[:TI])" transpose(optsol[:TI])
        # end

        optsol[:value] = problem[:costtransportation] * reduce(+, [min(sum(optsol[:TI][t, :]), 1) for t in 1:nbplannedtrucks], init=0.0) + 
        problem[:costextratruck] * reduce(+, [min(sum(optsol[:TI][t, :]), 1) for t in nbplannedtrucks+1:length(chosentrucks)], init=0.0) + 
        problem[:costinventory] * sum(problem[:IDL] - transpose(optsol[:TI]) * problem[:TDA][chosentrucks])
        
       
        clearnlines(1)
        printstyled("Solved with ", length(chosentrucks), " trucks\n", color=:light_blue)
        # If no unassigned item:
        if sum(optsol[:TI]) == nbitems
            # remove empty unused useless trucks
            for (i, t) in enumerate(chosentrucks)
                if sum(optsol[:TI][i, :]) == 0 && !(t in 1:nbplannedtrucks)
                    push!(deadtrucks, t)
                end
            end
            filter!(x -> !(x in deadtrucks), chosentrucks)
        end

        # improvement = optsol[:value] < oldvalue || optsol[:TI] != oldTI
        improvement = optsol[:value] < oldvalue + eps # - eps? TODO
        # @debug begin
        #     @debug "improvement" improvement
        #     @debug "optsol[:value]" optsol[:value]
        #     @debug "oldvalue" oldvalue
        #     @debug "optsol[:TI] != oldTI" optsol[:TI] != oldTI
        #     sleep(10)
        # end


    end
    printstyled("End of Column Generation\n", color=:light_green)
    printstyled("Value change: \n", color=:light_green)
    printstyled(optsol[:value] - oldvalue, "\n")
    printstyled("Number of unassigned items: \n", color=:light_green)
    printstyled(nbitems - sum(optsol[:TI]), "\n")
    printstyled("Solved with ", length(chosentrucks), " trucks ($nbitems items)\n", color=:light_blue)

    return optsol

end