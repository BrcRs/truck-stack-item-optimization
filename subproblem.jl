using JuMP
using MathOptInterface

using Base.Threads
# using Clp # Only for debug
# using CDD # Only debug
# include("instance_loader.jl")
include("matrix_ops.jl")
# include("linear_infeasibilities.jl")
# include("progress.jl")
"""
    Subproblem(t::Integer, problem::TSIProblem, submodel::Union{Model, Nothing}, valueTI::Matrix{Union{Missing, Bool}})

Store the data of a subproblem tied to truck `t`. Problem data is accessible in `problem`.
`valueTI` stores the truck-item solution matrix of the subproblem.

Data from `problem` can be accessed via the Subproblem. For instance, you can 
retrieve `problem[:TL]` via `subproblem[:TL]`.
"""
mutable struct Subproblem
    t::Integer
    problem::TSIProblem
    submodel::Union{Model, Nothing}
    # optimizer::Any # Optimizer
    # nbstacks::Integer
    valueTI::Matrix{Union{Missing, Bool}}
end
model(sub::Subproblem) = sub.submodel
optimizer(sub::Subproblem) = sub.optimizer
truck(sub::Subproblem) = sub.t
getproblem(sub::Subproblem) = sub.problem

function Base.getindex(sub::Subproblem, name::Symbol)
    obj_dict = object_dictionary(sub.problem)
    if !haskey(obj_dict, name)
        throw(KeyError(name))
    end
    return obj_dict[name]
end



function Base.haskey(subproblem::Subproblem, name::Symbol)
    return haskey(object_dictionary(subproblem.problem), name)
end
valueTI(subpb::Subproblem) = subpb.valueTI


function setvalueTI!(subpb::Subproblem, value::BitMatrix)
    subpb.valueTI .= value
end

"""
    upd_valueTI!(subpb::Subproblem)

Updates subpb.valueTI with `TI` of subpb.model. To do after solving the model.
"""
function upd_valueTI!(subpb::Subproblem)
    setvalueTI!(subpb, value.(model(subpb)[:TI]))
end

function suppr!(model::Model, key::Symbol)
    delete.(model, model[key])
    JuMP.unregister.(model, key)
end

## add_var! is unused
# function add_var!(model, key::Symbol, len::Integer; replace=false)
#     return add_var!(model, key::Symbol, (len); replace)
# end

# function add_var!(model, key::Symbol, lens::Tuple; replace=false)
#     if replace
#         delete.(model, model[key])
#         JuMP.unregister.(model, key)
#     end
#     if length(lens) == 1
#         model[key] = @variable(model, [1:lens[1]]; base_name=key)
#     elseif length(lens) == 2
#         model[key] = @variable(model, [1:lens[1], 1:lens[2]]; base_name=key)
#     elseif length(lens) == 3
#         model[key] = @variable(model, [1:lens[1], 1:lens[2], 1:lens[3]]; base_name=key)
#     else
#         error("length of $length(lens) for lens is not supported.")
#     end
    
#     return model[key]
# end

# function doif(cond, f1, f2, args...)
#     if cond
#         return f1(args...)
#     else
#         return f2(args...)
#     end
# end

"""
    buildTSImodel!(submodel::Model, problem::TSIProblem, t::Integer, chosentrucks::Vector{Integer}; replace=false)

Build a JuMP model to solve the truck stack item assignation subproblem for 
truck t to solve the larger problem the ROADEF challenge of 2022, only using `chosentrucks`.
Modify the model in place. If `replace=true`, will delete some constraints exclusive to `t` to replace them.
You want to keep all stack related constraints when the given submodel to update had a 
compatible truck with the new one.
"""
function buildTSImodel!(submodel::Model, problem::TSIProblem, t::Integer, chosentrucks::Vector{Integer}; replace=false)

    set_silent(submodel)
    # index in subproblem might differ from global indexing used in chosentrucks.
    # chosentrucks holds indices based on all planned trucks + extra trucks, while
    # here we have a mixture of non-contiguous planned trucks + extra trucks
    i_t = findfirst(==(t), chosentrucks)
    nbcandidateitems = sum(problem[:TR][t, :]) # TR is global matrix of all planned trucks + extra trucks, hence t is used
    # nbcandidateitems = max([sum(problem[:TR][t2, :]) for t2 in 1:nbtrucks]...)
    nbstacks =nbcandidateitems
    nbitems = problem[:nbitems]
    nbtrucks = problem[:nbtrucks]
    nbchosentrucks = length(chosentrucks)
    nbplannedtrucks = problem[:nbplannedtrucks]
    # nbplants = problem[:nbplants]
    nbplantdocks = problem[:nbplantdocks]
    nbsuppliers = problem[:nbsuppliers]
    nbsupplierdocks = problem[:nbsupplierdocks]
    nbextratrucks = length(filter(x -> !(x in 1:nbplannedtrucks), chosentrucks))
    ## Add variables
    # # @info "Creating variables..."
    # @info "Adding zetaT..."
    if !replace
        # zeta is used to determine if a truck is used or not
        @variable(submodel, zetaT[1:nbplannedtrucks] >= 0) 
        @variable(submodel, zetaE[1:nbextratrucks] >= 0)
        # TI: truck-item binary matrix
        @variable(submodel, TI[1:nbchosentrucks, 1:nbitems], lower_bound = 0, upper_bound = 1, Bin)
        # GI: binary vector of unassigned items
        # necessary relaxation for use of column generation
        @variable(submodel, GI[1:nbitems], lower_bound = 0, upper_bound = 1, Bin)
        # R: linearization variable
        @variable(submodel, R[1:nbitems, 1:nbsuppliers], lower_bound = 0, upper_bound = 1, Bin)
        @variable(submodel, Theta[1:nbitems, 1:nbsuppliers], lower_bound = 0, upper_bound = 1, Bin)
        # variable orientations of items, sometimes fixed when forced orientation is provided
        @variable(submodel, IOV[1:nbitems], lower_bound = 0, upper_bound = 1, Bin)
    end
    # Here we only replace some variables when the new truck is not of the same kind as the previous one
    # When we deal with an extra truck and the previous truck was of the same kind, some variables don't need to be replaced
    if !replace || !haskey(submodel, :S)
        # We don't add those variables when there are no candidate items for the truck
        if nbcandidateitems > 0
            # stackability code of each stack
            @variable(submodel, SS[1:nbstacks] >= 0)
            # plant of each stack 
            @variable(submodel, SP[1:nbstacks], lower_bound = 0)
            # plant dock of each stack
            @variable(submodel, SK[1:nbstacks, 1:nbsupplierdocks], lower_bound = 0, upper_bound = 1) # No need for Bin because it is constrained to be equal to their items' which are integer
            # Supplier (dock?) of each stack
            @variable(submodel, SU[1:nbstacks, 1:nbsuppliers], lower_bound = 0, upper_bound = 1)
            # orientation of each stack
            @variable(submodel, SO[1:nbstacks] >= 0)
            # X extremity of each stack
            @variable(submodel, SXe[1:nbstacks] >= 0)
            # X origin of each stack
            @variable(submodel, SXo[1:nbstacks] >= 0)
            # Y extremity of each stack
            @variable(submodel, SYe[1:nbstacks] >= 0)
            # etc
            @variable(submodel, SYo[1:nbstacks] >= 0)
            @variable(submodel, SZe[1:nbstacks] >= 0)
            # those are linearization variables
            @variable(submodel, betaM[1:convert(Int, nbstacks*(nbstacks+1)/2)] >= 0)
            @variable(submodel, betaP[1:convert(Int, nbstacks*(nbstacks+1)/2)] >= 0)
            # @variable(submodel, nu[1:nbstacks] >= 0)
            # @variable(submodel, tau[1:nbstacks] >= 0)
            # @variable(submodel, phi[1:nbstacks] >= 0)
            # Supplier dock of each stack
            @variable(submodel, SG[1:nbstacks, 1:nbplantdocks] >= 0, upper_bound = 1)
            # Length and width
            @variable(submodel, SL[1:nbstacks] >= 0)
            @variable(submodel, SW[1:nbstacks] >= 0)
            # stack-item assignment
            @variable(submodel, S[1:nbstacks, 1:nbcandidateitems], lower_bound = 0, upper_bound = 1, Bin)

            @variable(submodel, Z[1:nbstacks, 1:nbcandidateitems], lower_bound = 0, upper_bound = 1, Bin)

            ## many linearization variables
            @variable(submodel, mu[1:convert(Int, nbstacks*(nbstacks+1)/2)], lower_bound = 0, upper_bound = 1, Bin)
            # @variable(submodel, eta[1:convert(Int, nbstacks*(nbstacks+1)/2)], lower_bound = 0, upper_bound = 1, Bin)
            @variable(submodel, xi[1:convert(Int, nbstacks*(nbstacks+1)/2)], lower_bound = 0, upper_bound = 1, Bin)
            @variable(submodel, chi[1:convert(Int, nbstacks*(nbstacks+1)/2)], lower_bound = 0, upper_bound = 1, Bin)
            @variable(submodel, r[1:convert(Int, nbstacks*(nbstacks+1)/2)], lower_bound = 0, upper_bound = 1, Bin)
            @variable(submodel, sigma1[1:convert(Int, nbstacks*(nbstacks+1)/2)], lower_bound = 0, upper_bound = 1, Bin)
            @variable(submodel, sigma2[1:convert(Int, nbstacks*(nbstacks+1)/2)], lower_bound = 0, upper_bound = 1, Bin)
            @variable(submodel, sigma3[1:convert(Int, nbstacks*(nbstacks+1)/2)], lower_bound = 0, upper_bound = 1, Bin)
            @variable(submodel, Q[1:nbstacks, 1:nbsuppliers], lower_bound = 0, Int)
            @variable(submodel, V[1:nbstacks, 1:nbsupplierdocks], lower_bound = 0, Int)
            @variable(submodel, W[1:nbstacks, 1:nbplantdocks], lower_bound = 0, Int)
            @variable(submodel, Gl[1:nbstacks, 1:nbcandidateitems], lower_bound = 0, Int)
            @variable(submodel, Gr[1:nbstacks, 1:nbcandidateitems], lower_bound = 0, Int)
            @variable(submodel, DL[1:nbstacks, 1:nbcandidateitems], lower_bound = 0)
            @variable(submodel, DW[1:nbstacks, 1:nbcandidateitems], lower_bound = 0)
            @variable(submodel, lambda[1:convert(Int, nbstacks * (nbstacks+1)/2)], lower_bound = 0, Bin)
    
        end
    end
    ### begin

    ## Big M constants used in many linearizations
    MZ = max(problem[:IS]...) + 1.0
    MQ = max(problem[:IU]...) + 1.0
    MV = max(problem[:IK]...) + 1.0
    MW = max(problem[:IPD]...) + 1.0
    MG = 5.0
    MDL =  max(problem[:IL]...) + 1.0
    MDW =  max(problem[:IW]...) + 1.0
    MTL = max(problem[:TL]...) + 1.0
    MTW = max(problem[:TW]...) + 1.0
    MTE = max(skipmissing(problem[:TE])...)
    MTKE = max(problem[:TKE]...)
    MTGE = max(problem[:TGE]...)
    MTW = max(problem[:TW]...) + 1.0
    Mlambda = 2.0 * max(problem[:TL]...) + 1.0
    Mzeta = nbitems
    Mmu = 2.0
    # MS = 2.0

    # Those matrices help selecting every couple of stack
    # Xi1 is the first stack and Xi2 the second
    Xi1 = falses(convert(Int, nbstacks*(nbstacks+1)/2), nbstacks)
    Xi2 = falses(convert(Int, nbstacks*(nbstacks+1)/2), nbstacks)
    fillXi1!(Xi1)
    fillXi2!(Xi2)
    
    epsilon = 0.001

    ## Add constraints
    # When t changes and we reuse a model, some constraints need to be replaced
    if replace
        if haskey(JuMP.object_dictionary(submodel), :cZetaT2)
            suppr!(submodel, :cZetaT2)
        end
        if haskey(JuMP.object_dictionary(submodel), :cZetaT3)
            suppr!(submodel, :cZetaT3)
        end
        if haskey(JuMP.object_dictionary(submodel), :cZetaE2)
            suppr!(submodel, :cZetaE2)
        end
        if haskey(JuMP.object_dictionary(submodel), :cZetaE3)
            suppr!(submodel, :cZetaE3)
        end
    end
    if t <= nbplannedtrucks # TODO something is wrong here. Why the condition? Is not dependent on t.
        # sum(TI[t, :]) implies zeta[t] = 1. Else since we minimize, zeta = 0 
        @constraint(submodel, cZetaT2, submodel[:zetaT] * Mzeta .>= submodel[:TI][1:nbplannedtrucks, :] * vones(Int8, nbitems))
        @constraint(submodel, cZetaT3, -(1 .- submodel[:zetaT]) * Mzeta .+ 1 .<= submodel[:TI][1:nbplannedtrucks, :] * vones(Int8, nbitems))
    else
        @constraint(submodel, cZetaE2, submodel[:zetaE] * Mzeta .>= submodel[:TI][nbplannedtrucks+1:end, :] * vones(Int8, nbitems))
        @constraint(submodel, cZetaE3, -(1 .- submodel[:zetaE]) * Mzeta .+ 1 .<= submodel[:TI][nbplannedtrucks+1:end, :] * vones(Int8, nbitems))
    end
    
    if replace
        if haskey(submodel, :cTI_TR)
            suppr!(submodel, :cTI_TR)
        end
        if haskey(submodel, :cS_TI)
            suppr!(submodel, :cS_TI)
        end
        if haskey(submodel, :cSXe_ST_TL)
            suppr!(submodel, :cSXe_ST_TL)
        end
        if haskey(submodel, :cSYe_ST_TW)
            suppr!(submodel, :cSYe_ST_TW)
        end
        if haskey(submodel, :cSZe_ST_TH)
            suppr!(submodel, :cSZe_ST_TH)
        end
    end
    # Only add constraint for not all but only candidate items of truck t
    icandidates = findall((x) -> x == 1, problem[:TR][t, :])

    # Filter out non-candidate items from the solution
    @constraint(submodel, cTI_TR, submodel[:TI][i_t, :] .<= problem[:TR][t, :])
    # @debug begin
    #     @debug "i_t" i_t
    #     @debug "transpose(submodel[:S])" transpose(submodel[:S])
    #     @debug "vones(Int8, nbstacks)" vones(Int8, nbstacks)
    #     @debug "submodel[:TI][i_t, filter(x -> x in icandidates, 1:nbitems)]" submodel[:TI][i_t, filter(x -> x in icandidates, 1:nbitems)]

    # end
    if nbcandidateitems > 0
        # If an item is in a stack, the item is in the truck, and vice-versa
        @constraint(submodel, cS_TI, transpose(submodel[:S]) * vones(Int8, nbstacks) .== submodel[:TI][i_t, filter(x -> x in icandidates, 1:nbitems)])
        # stack dimensions don't exceed truck's dimensions
        @constraint(submodel, cSXe_ST_TL, submodel[:SXe] .<= problem[:TL][t])
        @constraint(submodel, cSYe_ST_TW, submodel[:SYe] .<= problem[:TW][t])
        @constraint(submodel, cSZe_ST_TH, submodel[:SZe] .<= problem[:TH][t])
    end

    if !replace
        # There is at most one truck per item, relaxed with GI # TODO shouldn't it be ==?
        @constraint(submodel, cTI_1_1, transpose(submodel[:TI])[filter(x -> x in icandidates, 1:nbitems), :] * vones(Int8, nbchosentrucks) + submodel[:GI][filter(x -> x in icandidates, 1:nbitems)] .<= vones(Int8, size(transpose(submodel[:TI])[filter(x -> x in icandidates, 1:nbitems), :], 1)))
        if nbcandidateitems > 0
            # All items from a stack have same stackability code
            @constraint(submodel, cZ_S_MZ, submodel[:Z] .<= submodel[:S] * MZ)
            @constraint(submodel, cS_IS_Z, submodel[:S] * problem[:IS][filter(x -> x in icandidates, 1:nbitems)] .== submodel[:Z] * vones(Int8, size(submodel[:Z])[1]))
            @constraint(submodel, cZ_SS_Sleft, -MZ * (-submodel[:S] .+ 1) .<= submodel[:Z] .- hcat([submodel[:SS] for i in 1:size(submodel[:Z])[2]]...))
            @constraint(submodel, cZ_SS_Sright, submodel[:Z] .- hcat([submodel[:SS] for i in 1:size(submodel[:Z])[2]]...) .<= MZ * (-submodel[:S] .+ 1))
            
            # All items from a stack have same supplier
            @constraint(submodel, cQ_SU, submodel[:Q] .<= submodel[:SU] * MQ)
            @constraint(submodel, cS_IU_Q, submodel[:S] * problem[:IU][filter(x -> x in icandidates, 1:nbitems)] .== submodel[:Q])
            @constraint(submodel, cQ_SU_Sleft, -MQ * (1 .- submodel[:SU]) .<= submodel[:Q] .- hcat([submodel[:S] * vones(Int8, size(submodel[:S], 2)) for i in 1:size(submodel[:Q], 2)]...))
            @constraint(submodel, cQ_SU_Sright, submodel[:Q] .- hcat([submodel[:S] * vones(Int8, size(submodel[:S], 2)) for i in 1:size(submodel[:Q], 2)]...) .<= MQ * (1 .- submodel[:SU]))

            # All items from a stack must have the same supplier dock
            @constraint(submodel, cV_SK, submodel[:V] .<= submodel[:SK] * MV)
            @constraint(submodel, cS_IK_V, submodel[:S] * problem[:IK][filter(x -> x in icandidates, 1:nbitems)] .== submodel[:V])
            @constraint(submodel, cV_SK_Sleft, -MV * (-submodel[:SK] .+ 1) .<= submodel[:V] .- hcat([submodel[:S] * vones(Int8, size(submodel[:S], 2)) for i in 1:size(submodel[:V], 2)]...))
            @constraint(submodel, cV_SK_Sright, submodel[:V] .- hcat([submodel[:S] * vones(Int8, size(submodel[:S], 2)) for i in 1:size(submodel[:V], 2)]...) .<= MV * (-submodel[:SK] .+ 1))

            # Same with plant dock
            @constraint(submodel, cW_SG, submodel[:W] .<= submodel[:SG] * MW)
            @constraint(submodel, cS_IPD_W, submodel[:S] * problem[:IPD][filter(x -> x in icandidates, 1:nbitems)] .== submodel[:W])
            @constraint(submodel, cW_SPD_Sleft, -MW * (-submodel[:SG] .+ 1) .<= submodel[:W] .- hcat([submodel[:S] * vones(Int8, size(submodel[:S], 2)) for i in 1:size(submodel[:W], 2)]...))
            @constraint(submodel, cW_SPD_Sright, submodel[:W] .- hcat([submodel[:S] * vones(Int8, size(submodel[:S], 2)) for i in 1:size(submodel[:W], 2)]...) .<= MW * (-submodel[:SG] .+ 1))

            # Same with orientation. More complicated because orientation is a variable.
            @constraint(submodel, cGl_S, submodel[:Gl] .<= submodel[:S] * MG)
            @constraint(submodel, cGr_S, submodel[:Gr] .<= submodel[:S] * MG)
            @constraint(submodel, cGl_Gr, submodel[:Gl] .== submodel[:Gr])
            @constraint(submodel, cGr_SO_Sleft, -MG * (-submodel[:S] .+ 1) .<= submodel[:Gr] .- hcat([submodel[:SO] for i in 1:nbcandidateitems]...))
            @constraint(submodel, cGr_SO_Sright, submodel[:Gr] .- hcat([submodel[:SO] for i in 1:nbcandidateitems]...) .<= MG * (-submodel[:S] .+ 1))
            @constraint(submodel, cGl_IOV_Sleft, -MG * (-submodel[:S] .+ 1) .<= submodel[:Gl] .- vcat([transpose(submodel[:IOV][filter(x -> x in icandidates, 1:nbitems)]) for i in 1:nbstacks]...))
            @constraint(submodel, cGl_IOV_Sright, submodel[:Gl] .- vcat([transpose(submodel[:IOV][filter(x -> x in icandidates, 1:nbitems)]) for i in 1:nbstacks]...) .<= MG * (-submodel[:S] .+ 1))

            # All items of a stack have same length and width
            @constraint(submodel, cDL_S, submodel[:DL] .<= submodel[:S] * MDL)
            @constraint(submodel, cDL_SLleft, -MDL * (-submodel[:S] .+ 1) .<= submodel[:DL] - hcat([submodel[:SL] for i in 1:size(submodel[:DL])[2]]...))
            @constraint(submodel, cDL_SLright, submodel[:DL] - hcat([submodel[:SL] for i in 1:size(submodel[:DL])[2]]...) .<= MDL * (-submodel[:S] .+ 1))
            @constraint(submodel, cDL_S_IL, submodel[:DL] * vones(Int8, size(submodel[:S], 1)) .== submodel[:S] * problem[:IL][filter(x -> x in icandidates, 1:nbitems)])
            @constraint(submodel, cDW_S, submodel[:DW] .<= submodel[:S] * MDW)
            @constraint(submodel, cDW_SWleft, -MDW * (-submodel[:S] .+ 1) .<= submodel[:DW] - hcat([submodel[:SW] for i in 1:size(submodel[:DW])[2]]...))
            @constraint(submodel, cDW_SWright, submodel[:DW] - hcat([submodel[:SW] for i in 1:size(submodel[:DW])[2]]...) .<= MDW * (-submodel[:S] .+ 1))
            @constraint(submodel, cDW_S_IW, submodel[:DW] * vones(Int8, size(submodel[:S], 1)) .== submodel[:S] * problem[:IW][filter(x -> x in icandidates, 1:nbitems)])

            # Define coordinates of each stack based on dimension and orientation
            @constraint(submodel, cSXe_SXo_SL_SOleft, submodel[:SXe] - submodel[:SXo] - submodel[:SL] .<= submodel[:SO] * MTL)
            @constraint(submodel, cSXe_SXo_SL_SOright, submodel[:SXe] - submodel[:SXo] - submodel[:SL] .>= -submodel[:SO] * MTL)
            @constraint(submodel, cSYe_SYo_SW_SOleft, submodel[:SYe] - submodel[:SYo] - submodel[:SW] .<= submodel[:SO] * MTW)
            @constraint(submodel, cSYe_SYo_SW_SOright, submodel[:SYe] - submodel[:SYo] - submodel[:SW] .>= -submodel[:SO] * MTW)
            @constraint(submodel, cSXe_SXo_SW_SOleft, submodel[:SXe] - submodel[:SXo] - submodel[:SW] .<= (-submodel[:SO] .+ 1) * MTW)
            @constraint(submodel, cSXe_SXo_SW_SOright, submodel[:SXe] - submodel[:SXo] - submodel[:SW] .>= -(-submodel[:SO] .+ 1) * MTW)
            @constraint(submodel, cSYe_SYo_SL_SOleft, submodel[:SYe] - submodel[:SYo] - submodel[:SL] .<= (-submodel[:SO] .+ 1) * MTL)
            @constraint(submodel, cSYe_SYo_SL_SOright, submodel[:SYe] - submodel[:SYo] - submodel[:SL] .>= -(-submodel[:SO] .+ 1) * MTL)

            # Compute height of stack
            @constraint(submodel, cSZe_S_IH, submodel[:SZe] .== submodel[:S]* problem[:IH][filter(x -> x in icandidates, 1:nbitems)])

            # Set an order between all stacks, making placement constraints simpler
            @constraint(submodel, cSXo_SXo, submodel[:SXo][1:end-1] .<= submodel[:SXo][2:end])

            # Determine which of two stacks precedes the other
            # If they overlap on X, they won't overlap on Y due to betaM which determines mu which determines Y constraints further down
            @constraint(submodel, cXi2SXo_Xi1SXe_betaM_betaP, (Xi2 * submodel[:SXo]) - (Xi1 * submodel[:SXe]) - submodel[:betaM] + submodel[:betaP] .== -epsilon * vones(Float64, size(Xi1, 1)))

            @constraint(submodel, cbetaM_lambda, submodel[:betaM] .<= submodel[:lambda] .* Mlambda)
            @constraint(submodel, cbetaP_lambda, submodel[:betaP] .<= (-submodel[:lambda] .+ 1)*Mlambda)
            @constraint(submodel, cmu_betaM, (-submodel[:mu] .+ 1) .<= submodel[:betaM] * Mmu)

            # Stacks overlap on Y axis if and only if these aren't overlapping on X axis (variable mu)
            # xi implies only one of two constraints must be satisfied
            @constraint(submodel, cXi1SYe_Xi2SYo, Xi1 * submodel[:SYe] .<= Xi2 * submodel[:SYo] + submodel[:xi] * MTW + (-submodel[:mu] .+ 1) * MTW)
            @constraint(submodel, cXi2SYe_Xi1SYo, Xi2 * submodel[:SYe] .<= Xi1 * submodel[:SYo] + (-submodel[:xi] .+ 1) * MTW + (-submodel[:mu] .+ 1) * MTW)

            ## Determine pick-up/loading orders for supplier, supplier dock, plant dock
            notmissingTE = filter(x -> !ismissing(problem[:TE][t, x]), 1:nbsuppliers)
            @constraint(submodel, cXi1SU_Xi2SU, Xi1 * submodel[:SU][:, notmissingTE] * problem[:TE][t, notmissingTE] .<= Xi2 * submodel[:SU][:, notmissingTE] * problem[:TE][t, notmissingTE])
            @constraint(submodel, cXi1SU_Xi2SU_chi, Xi1*submodel[:SU] - Xi2*submodel[:SU] .>= submodel[:chi] * epsilon - submodel[:r]*MTE - (-submodel[:sigma1] .+ 1) * MTE)
            @constraint(submodel, cXi2SU_Xi1SU_chi, Xi2*submodel[:SU] - Xi1*submodel[:SU] .>= (-submodel[:chi] .+ 1) * epsilon - submodel[:r]*MTE - (-submodel[:sigma1] .+ 1) * MTE)

            notmissingTKE = filter(x -> !ismissing(problem[:TKE][t, x]), 1:nbsupplierdocks)
            @constraint(submodel, cXi2SK_Xi1SK, Xi2*submodel[:SK][:, notmissingTKE]*problem[:TKE][t, notmissingTKE] .>= Xi1*submodel[:SK][:, notmissingTKE]*problem[:TKE][t, notmissingTKE] - (-submodel[:r] .+ 1) * MTKE)
            @constraint(submodel, cXi1SK_Xi2SK_chi, Xi1*submodel[:SK][:, notmissingTKE]*problem[:TKE][t, notmissingTKE] - Xi2*submodel[:SK][:, notmissingTKE]*problem[:TKE][t, notmissingTKE] .>= submodel[:chi]*epsilon - (-submodel[:sigma2] .+ 1)*MTKE)
            @constraint(submodel, cXi2SK_Xi1SK_chi, Xi2*submodel[:SK][:, notmissingTKE]*problem[:TKE][t, notmissingTKE] - Xi1*submodel[:SK][:, notmissingTKE]*problem[:TKE][t, notmissingTKE] .>= (-submodel[:chi] .+ 1)*epsilon - (-submodel[:sigma2] .+ 1)*MTKE)

            notmissingTGE = filter(x -> !ismissing(problem[:TGE][t, x]), 1:nbplantdocks)
            @constraint(submodel, cXi2SG_Xi1SG, Xi2*submodel[:SG][:, notmissingTGE]*problem[:TGE][t, notmissingTGE] .>= Xi1*submodel[:SG][:, notmissingTGE]*problem[:TGE][t, notmissingTGE] - (-submodel[:sigma3] .+ 1) * MTGE)

            @constraint(submodel, csigma1_sigma2_sigma3, submodel[:sigma1] + submodel[:sigma2] + submodel[:sigma3] .>= 1)
        end
    end

    return submodel
end

function Subproblem(t, problem::TSIProblem, chosentrucks)
    submodel = buildTSImodel!(Model(optimizer(problem)), problem, t, chosentrucks)
    subproblem = Subproblem(t, problem, submodel, Matrix{Union{Missing, Bool}}(missing, problem[:nbtrucks], problem[:nbtrucks]))
    return subproblem
end

function setmodel!(subproblem::Subproblem, model)
    subproblem.submodel = model
end