using JuMP
using MathOptInterface

using Base.Threads
# using Clp # Only for debug
# using CDD # Only debug
include("instance_loader.jl")
include("matrix_ops.jl")
# include("linear_infeasibilities.jl")
# include("progress.jl")

"""
    TSIProblem(opt, obj_dict::Dict{Symbol,Any})

Problem data of a truck stack item affectation problem.
`opt` is an optimizer. `obj_dict` holds all problem data.

You can directly access the dictionary with problem[:something].
"""
mutable struct TSIProblem
    opt::Any # Optimizer
    obj_dict::Dict{Symbol,Any}
end

object_dictionary(pb::TSIProblem) = pb.obj_dict

optimizer(pb::TSIProblem) = pb.opt

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

