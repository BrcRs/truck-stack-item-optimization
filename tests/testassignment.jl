using Test

include("../src/assignment.jl")
include("../src/instance_loader.jl")

@testset "solve_tsi" begin
    instancepath = "./instances/AS/"
    item_productcodes, truckdict, supplierdict, supplierdockdict, plantdict, 
    plantdockdict, nbplannedtrucks, nbitems, nbsuppliers, nbsupplierdocks, nbplants, nbplantdocks,
    TE_P, TL_P, TW_P, TH_P, TKE_P, TGE_P, TDA_P, TDE_P, TU_P, TP_P, TK_P, TG_P, TR_P,
    IU, IP, IK, IPD, IS, _IO, IL, IW, IH, IDL, IDE, stackabilitycodedict, nbtrucks, TE, TL, TW, TH,
    TKE, TGE, TDA, TDE, TU, TP, TK, TG, TR, TID, reverse_truckdict, truckindices,
    costinventory, costtransportation, costextratruck, timelimit = loadinstance(instancepath; onlyplanned=true)
    
    error("TEM = ???") # TEM is max stack density
    error("TMm = ???")

    display(size(TR_P))

    t_trucks = []
    for t in 1:nbplannedtrucks
        push!(
            t_trucks, 
            Pair(
                t, 
                Truck(
                    reverse_truckdict[t],
                    Dim(TL_P[t], TW_P[t]),
                    TH_P[t],
                    TEM_P[t],
                    Dict{Any, Any}(),
                    TMm_P[t],
                    Dict(s => TE_P[t, supplierdict[s]] for s in keys(supplierdict) if TU_P[t, supplierdict[s]]),
                    Dict(s => Dict(sd => TKE_P[t, supplierdockdict[sd]] for sd in keys(supplierdockdict[s])) 
                        for s in keys(supplierdict) if TU_P[t, supplierdict[s]]),
                    Dict(p => Dict(pd => TGE_P[t, plantdockdict[pd]] for pd in keys(plantdockdict[p])) 
                        for p in keys(plantdict) if TP_P[t, plantdict[p]]),
                    costextratruck,
                    error("Continue here")

                )
            )
        )
    end
    # TODO set max_stack_weights

    # solve_tsi(t_trucks, i_items, TR)
    
end