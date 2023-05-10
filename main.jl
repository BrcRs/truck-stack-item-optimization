
using LinearAlgebra
using JuMP
using Cbc, GLPK, Clp
# using CSV

# using FilePaths
using Logging
# using Profile

using ArgParse


include("instance_loader.jl")
include("model.jl")
include("column_generation.jl")
include("tsiproblem.jl")
include("uzawa.jl")

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
        "--debug", "-d"
            help= "activate debug loggings"
            action = :store_true
        "--all", "-a"
            help = "shows unlimited debug loggings"
            action = :store_true
        "--optimizer", "-o"
            help = "chosen optimizer"
            default = "Cbc"
            arg_type = String
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
    optimizers = Dict{String, Any}("Cbc" => Cbc, "GLPK" => GLPK)
    if parsed_args["debug"]
        logger = ConsoleLogger(stdout, Logging.Debug, show_limited=!parsed_args["all"])

        global_logger(logger)
    end
    # instancePath = "Instances/AS/"
    instancePath = parsed_args["instancePath"]
    opt = optimizers[parsed_args["optimizer"]]
    problem = TSIProblem(opt.Optimizer, instancePath)
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
    @time begin 
        optsol = columngeneration(solve_uzawa!, problem, 10, 10, 1)
    end
    display(optsol[:TI])
    selectedtrucks = filter(x -> sum(optsol[:TI][x, :]) >= 1, 1:size(optsol[:TI], 1))
    printstyled("Selected trucks:\n", color=:green)
    display([t <= plannedtrucks ? string("P", t) : string("E", t), selectedtrucks])
    for t in selectedtrucks
        print(t, ":")
        display(optsol[:variables][t])
    end
end
    main()
