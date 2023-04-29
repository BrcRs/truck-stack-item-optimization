using JuMP

include("instance_loader.jl")


function main()
    instancepath = "Instances/CA/"

    item_productcodes, truckdict, supplierdict, supplierdockdict, plantdict, 
    plantdockdict, nbplannedtrucks, nbitems, nbsuppliers, nbsupplierdocks, nbplants, nbplantdocks,
    TE_P, TL_P, TW_P, TH_P, TKE_P, TGE_P, TDA_P, TU_P, TP_P, TK_P, TG_P, TR_P,
    IU, IP, IK, IPD, IS, _IO, IL, IW, IH, IDL, IDE, stackabilitycodedict, nbtrucks, TE, TL, TW, TH,
    TKE, TGE, TDA, TU, TP, TK, TG, TR, TID, reverse_truckdict, 
    costinventory, costtransportation, costextratruck, timelimit = loadinstance(instancepath)

    # Put each item in a single empty matching truck as to minimize overall cost
    # this constitutes the root node of the search tree.

    # root = Set{??}()

    # subpb[:costtransportation] * submodel[:zetaT] + 
    # subpb[:costextratruck] * submodel[:zetaE] + 
    # subpb[:costinventory] * (subpb[:IDL] - transpose(submodel[:TI][t]) * subpb[:TDA]) + 
    # kappa * (submodel[:TI] - TIbar))

    root = [Dict{String, Set{Any}}(), 0]

    bannedt = Set{String}()

    emptytrucks = Set{String}([string(t, "P") for t in 1:nbplannedtrucks])

    extratruckcount = Dict(string(t, "E") => 0 for t in 1:nbplannedtrucks)

    @debug "nbitems" nbitems
    @debug "nbitems^2" nbitems * nbitems
    @debug "2nbitems - 1" 2*nbitems - 1

    for i in 1:nbitems
        cost = 0
        t = ""
        if !isempty(filter(x -> TR_P[parse(Int64, x[1:end-1]), i] == 1, emptytrucks))
            t = argmin(t -> costtransportation + costinventory * (IDL[i] - TDA[parse(Int64, t[1:end-1])]), 
            filter(x -> TR_P[parse(Int64, x[1:end-1]), i] == 1, emptytrucks))
            cost = costtransportation + costinventory * (IDL[i] - TDA[parse(Int64, t[1:end-1])])
            pop!(emptytrucks, t)
            ## ...
        else
            t = argmin(t -> costextratruck + costinventory * (IDL[i] - TDA[parse(Int64, t[1:end-1])]), 
            [string(t, "E") for t in filter(x -> TR_P[x, i] == 1, 1:nbplannedtrucks)])
            cost = costextratruck + costinventory * (IDL[i] - TDA[parse(Int64, t[1:end-1])])
            # @debug "bannedt" bannedt
            extratruckcount[t] += 1
            t = string(t, extratruckcount[t])
            # @debug "t" t
            # sleep(10) # @debug
            ## ...
        end
        push!(bannedt, t)
        if !(t in keys(root[1]))
            root[1][t] = Set{Any}()
        end
        push!(root[1][t], i)
        root[2] += cost

    end

    # show(root)

    # The number of neighbors for a given level of n trucks is given by n(n-1)/2
    # The number of level can be approximated by the number of nodes of a complete binary tree which is 2nbitems - 1
    # We thus have at most (nbitems^3 - nbitems^2)/2 iterations

    # For each neighbor of the current node, a neighbor consisting in moving an 
    # item from a truck to another:
    # Better: a neighbor is two trucks fusioning:

        # Estimate feasibility

        # If not memoized, Estimate lower and upper bound on cost
            # The upper bound can't be worse than its predecessor's
    # save best lower bound
    # delete each solution with a lower bound higher than the best (lowest) upper bound
    # add neighbors to frontier

    # Take the node from the frontier with best estimated upper bound cost as current node
    # Repeat

    ## Good upper bound estimator:
    # Price associated with arrival, truck cost will be exact

    ## Feasibility estimator:
    # If sum((S * IL)[s, :] <= TL: FEASIBLE

    # Put a price on feasibility?
    # Get upper bound by adding a major simplifying constraint
    # E.g. All stacks have same dimensions which is max
    # Get lower bound by relaxing a subtle complexifying constraint
    # E.g. Relax loading orders (might not be enough)
    # E.g. All stacks have same dimensions which is min

end

main()