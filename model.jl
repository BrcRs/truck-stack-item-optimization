using JuMP
using MathOptInterface

using Base.Threads
# using Clp # Only for debug
# using CDD # Only debug
include("instance_loader.jl")
include("matrix_ops.jl")
include("linear_infeasibilities.jl")
include("progress.jl")
# TODO Refactor
