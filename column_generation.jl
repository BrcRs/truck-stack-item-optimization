using LinearAlgebra

function columngeneration(solvefun!, problem, args...)

    nbplannedtrucks = problem[:nbplannedtrucks]
    truckindices = problem[:truckindices]
    deadtrucks = []
    chosentrucks = [truckindices[t][j] for t in 1:nbplannedtrucks for j in 1:2]
    firstextratrucks = [truckindices[t][1] for t in 1:nbplannedtrucks]
    optsol = solvefun!(problem, args..., chosentrucks)

    improvement = true

    oldTI = copy(optsol[:TI])
    oldvalue = optsol[:value] * 2


    while transpose(optsol[:TI]) * vones(Int8, nbtrucks) != vones(Int8, nbitems) || improvement

        # Upd chosen trucks
        # Add extra trucks for each filled planned trucks
        for t in chosentrucks
            if sum(optsol[:TI][t, :]) > 0
                if t in 1:nbplannedtrucks
                    if !(truckindices[t][2] in chosentrucks) && !(truckindices[t][2] in deadtrucks)
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


        optsol = solvefun!(problem, args..., chosentrucks)
        
        # If no unaffected item:
        if sum(optsol[:GI]) == 0
            # remove empty unused useless trucks
            for t in chosentrucks
                if sum(optsol[:TI][t, :]) == 0 && !(t in 1:nbplannedtrucks)
                    push!(deadtrucks, t)
                end
            end
            filter!(x -> x in deadtrucks, chosentrucks)
        end

        improvement = optsol[:value] < oldvalue || optsol[:TI] != oldTI
        oldTI = copy(optsol[:TI])
        oldvalue = optsol[:value]
    end

    return optsol

end