using CSV
using Tables

include("item.jl")

# need a function which takes a truck, and a solution of type Dict{Integer, ItemizedStack}
# and prints to files output_items.csv, output_stacks.csv and output_trucks.csv


function solution_to_csv(truck, solution, directory; append=true)

    truck_mat = Dict{String, Any}()
    stack_mat = Dict{String, Vector{Any}}()
    item_mat = Dict{String, Vector{Any}}()

    # truck: columns needed
    # - Id truck
    truck_mat["Id truck"] = get_id(truck)
    # - Loaded length
    truck_mat["Loaded length"] = get_loaded_length([p[2] for p in solution])
    # - Weight of loaded items
    truck_mat["Weight of loaded items"] = get_loaded_weight([p[2] for p in solution])
    # - Volume of loaded items
    truck_mat["Volume of loaded items"] = get_loaded_volume([p[2] for p in solution])
    
    tm_t, ej_e, ej_r, em_h, em_r, em_m = 
        dist_stacks_to_trailer([p[2] for p in solution], truck)
    # - EMm
    truck_mat["EMm"] = em_m
    # - EMr
    truck_mat["EMr"] = em_r
    # then four unnamed columns (EMr?)
    # TODO

    # stacks: columns needed
    stack_cols = [
        "Id truck", "Id stack", "Stack code", "X origin", "Y origin", 
        "Z origin", "X extremity", "Y extremity", "Z extremity"
    ]
    for col in stack_cols
        stack_mat[col] = Vector{Any}()
    end
    item_cols = [
        "Id ident", "Id truck", "Id stack", "Item code", "X origin", "Y origin", 
        "Z origin", "X extremity", "Y extremity", "Z extremity"
    ]
    for col in item_cols
        item_mat[col] = Vector{Any}()
    end
    for (i, stack) in solution
        # - Id truck
        push!(stack_mat["Id truck"], get_id(truck))
        # - Id stack
        if get_id(stack) == ""
            error("Stacks must have unique ids.")
        end
        push!(stack_mat["Id stack"], get_id(stack))
        # - Stack code
        push!(stack_mat["Stack code"], get_minmax_stackability(stack))
        # - X origin
        push!(stack_mat["X origin"], get_pos(stack).x)
        # - Y origin
        push!(stack_mat["Y origin"], get_pos(stack).y)
        # - Z origin
        push!(stack_mat["Z origin"], 0) # TODO ? It is not always 0?
        # - X extremity
        push!(stack_mat["X extremity"], get_pos(stack).x + get_dim(stack).le)
        # - Y extremity
        push!(stack_mat["Y extremity"], get_pos(stack).y + get_dim(stack).wi)
        # - Z extremity
        push!(stack_mat["Z extremity"], get_height(stack))

        for (j, item) in enumerate(get_items(stack))
            # items: columns needed
            # - Item ident 
            push!(item_mat["Id ident"], get_id(item))
            # - Id truck 
            push!(item_mat["Id truck"], get_id(truck))
            # - Id stack 
            push!(item_mat["Id stack"], get_id(stack))
            # - Item code 
            push!(item_mat["Item code"], get_stackability_code(item))
            # - X origin 
            push!(item_mat["X origin"], get_pos(stack).x)
            # - Y origin 
            push!(item_mat["Y origin"], get_pos(stack).y)
            # - Z origin 
            push!(item_mat["Z origin"], (j-1) * get_height(item))
            # - X extremity 
            push!(item_mat["X extremity"], get_pos(stack).x + get_dim(stack).le)
            # - Y extremity 
            push!(item_mat["Y extremity"], get_pos(stack).y + get_dim(stack).wi)
            # - Z extremity
            push!(item_mat["Z extremity"], j * get_height(item))
        end
    end



    # TODO convert to table
    truck_tab = Tables.dictrowtable(truck_mat)
    stack_tab = Tables.dictrowtable(stack_mat)
    item_tab = Tables.dictrowtable(item_mat)
    
    # display(truck_mat)
    # readline()
    # display(stack_mat)
    # readline()
    # display(item_mat)
    # readline()

    touch(joinpath(directory, "output_trucks.csv"))
    touch(joinpath(directory, "output_stacks.csv"))
    touch(joinpath(directory, "output_items.csv"))

    # TODO then write to directory
    # output_items.csv, output_stacks.csv and output_trucks.csv
    CSV.write(joinpath(directory, "output_trucks.csv"), truck_tab, append=append, writeheader=!append)
    CSV.write(joinpath(directory, "output_stacks.csv"), stack_tab, append=append, writeheader=!append)
    CSV.write(joinpath(directory, "output_items.csv"), item_tab, append=append, writeheader=!append)
    error("Columns are not in order")
    error("Truck file is messed up")
end


# # for testing purposes, I need a functoin which takes output files and creates
# # a Truck and a solution::Vector{Pair{Integer, ItemizedStack}}
# Impossible: solution files miss data

# function csv_to_solution(directory)
#     # CSV.File(input_itemsfile, normalizenames=true, delim=';', decimal=',', stripwhitespace=true, types=String)

#     CSV.File(joinpath(directory, "output_trucks.csv"))
#     CSV.File(joinpath(directory, "output_stacks.csv"))
#     CSV.File(joinpath(directory, "output_items.csv"))

# end