
using LinearAlgebra
using JuMP
using Cbc
# using CSV

# using FilePaths
# using Logging
# using Profile

using ArgParse


include("instance_loader.jl")
include("model.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        ## Arg Parse examples
        # "--opt1"
        #     help = "an option with an argument"
        # "--opt2", "-o"
        #     help = "another option with an argument"
        #     arg_type = Int
        #     default = 0
        # "--flag1"
        #     help = "an option without argument, i.e. a flag"
        #     action = :store_true
        "instancePath"
            help = "Path to the folder of the instance to solve."
            required = true
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end
    
    # instancePath = "Instances/AS/"
    instancePath = parsed_args["instancePath"]

    problem = TSIProblem(Cbc.Optimizer, instancePath)
    TID = problem[:TID]
    TR_P = problem[:TR_P]
    @debug begin
        open("tmp2.txt", "w") do io
            show(io, "text/plain", [string(problem[:reverse_truckdict][i], " : ", sum(TR_P[i,:])) for i in 1:size(TR_P)[1]])
        end

        
        open("tmp.txt", "w") do io
            show(io, "text/plain", TID)
        end
    end
    if sum(sum(TR_P[:, j]) == 0 ? 1 : 0 for j in 1:size(TR_P)[2]) > 0
        @warn string("Nb of items with no truck available according to TR: ", sum(sum(TR_P[:, j]) == 0 ? 1 : 0 for j in 1:size(TR_P)[2]))
    end

    @info "Solving problem..."
    solve_uzawa!(problem, 1, 1, 5)
end
    @allocated main()
