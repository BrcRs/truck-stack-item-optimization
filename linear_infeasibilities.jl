using JuMP
using MathOptInterface

const MOI = MathOptInterface

using Clp # Only for debug
using CDDLib # Only debug


function find_problematic_constraint!(originalmodel)
    @info "Let's find infeasible constraints!"
    @warn "The model's constraints are going to be changed"
    model, referencemap = copy_model(originalmodel)
    # @debug "referencemap" referencemap
    set_optimizer(model, Clp.Optimizer)
    set_silent(model)
    set_optimizer_attribute(model, "LogLevel", 0)
    set_optimizer_attribute(model, "MaximumIterations", 3)
    set_optimizer_attribute(model, "InfeasibleReturn", 1)
    set_time_limit_sec(model, 60)
    infeasible_constrs = []
    constraint_types = list_of_constraint_types(originalmodel)
    for t in constraint_types
        if t in list_of_constraint_types(originalmodel)
            for c in all_constraints(originalmodel, t[1], t[2])
                modelcopy, refmapcopy = copy_model(model)
                # f = MOI.get(originalmodel, MOI.ConstraintFunction(), c)
                
                # name = MOI.get(originalmodel, MOI.ConstraintName(), c)
                # @debug "f" f
                # @debug "typeof(f)" typeof(f)
                # @debug "name" name
                # error("Stoooop")
                delete(modelcopy, refmapcopy[referencemap[c]])
                # JuMP.unregister(model, c)
                try
                    optimize!(modelcopy)
                catch LoadError
                end
                if termination_status(modelcopy) != INFEASIBLE
                    @info "Problematic constraint found" c
                    push!(infeasible_constrs, c)
                end
                # model[Symbol(name)] = @constraint(model, f)
            end
        end
    end
    # @info "No problematic constraint found"
    return infeasible_constrs
end