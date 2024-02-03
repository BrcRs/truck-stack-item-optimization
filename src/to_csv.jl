using CSV
using Tables
using OrderedCollections

include("item.jl")

# need a function which takes a truck, and a solution of type Dict{Integer, ItemizedStack}
# and prints to files output_items.csv, output_stacks.csv and output_trucks.csv


function solution_to_csv(truck, solution, directory; append=true)

    truck_mat = OrderedDict{String, Vector{Any}}()
    stack_mat = OrderedDict{String, Vector{Any}}()
    item_mat = OrderedDict{String, Vector{Any}}()

    truck_cols = [
        "Id truck", "Loaded length", "Weight of loaded items", "Volume of loaded items", "EMm", 
        "EMr"
    ]
    for col in truck_cols
        truck_mat[col] = Vector{Any}()
    end

    # truck: columns needed
    # - Id truck
    push!(truck_mat["Id truck"], get_id(truck))
    # - Loaded length
    push!(truck_mat["Loaded length"], get_loaded_length([p[2] for p in solution]))
    # - Weight of loaded items
    push!(truck_mat["Weight of loaded items"], get_loaded_weight([p[2] for p in solution]))
    # - Volume of loaded items
    push!(truck_mat["Volume of loaded items"], get_loaded_volume([p[2] for p in solution]))
    
    tm_t, ej_e, ej_r, em_h, em_r, em_m = 
        dist_stacks_to_trailer([p[2] for p in solution], truck)
    # - EMm
    push!(truck_mat["EMm"], em_m)
    # - EMr
    push!(truck_mat["EMr"], em_r)
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
    # error("Columns are not in order")
    # error("Truck file is messed up")
end

function write_input(trucks, items, directory)

    truck_mat = OrderedDict{String, Vector{Any}}()
    item_mat = OrderedDict{String, Vector{Any}}()


    truck_cols = [
        "Item ident", "Supplier code", "Supplier dock", "Plant code", "Plant dock", 
        "Product code", "Package code", "Number of items", "Length", "Width", "Height", "Weight",
        "Nesting height", "Stackability code", "Forced orientation", "Earliest arrival time",
        "Latest arrival time", "Inventory cost", "Max stackability"
    ]
    for col in truck_cols
        truck_mat[col] = Vector{Any}()
    end

    item_cols = [
        "Supplier code", "Supplier loading order", "Supplier dock", "Supplier dock loading order",
        "Plant code", "Plant dock", "Plant dock loading order", "Product code", "Arrival time",
        "Id truck", "Length", "Width", "Height", "Max weight", "Stack with multiple docks",
        "Max density", "Max weight on the bottom item in stacks", "Cost", "EMmm", "EMmr", "CM"
        "CJfm", "CJfc", "CJfh", "EM", "EJhr", "EJcr", "EJeh"
    ]
    for col in item_cols
        item_mat[col] = Vector{Any}()
    end

    # write item infos
    # Item ident, Supplier code, Supplier dock, Plant code, Plant dock, 
    # Product code, Package code, Number of items, Length, Width, Height, Weight,
    # Nesting height, Stackability code, Forced orientation, Earliest arrival time,
    # Latest arrival time, Inventory cost, Max stackability

    for item in items
        push!(item_mat["Item ident"], get_id(item))
        push!(item_mat["Supplier code"], get_supplier(item))
        push!(item_mat["Supplier dock"], get_supplier_dock(item))
        push!(item_mat["Plant code"], get_plant(item))
        push!(item_mat["Plant dock"], get_plant_dock(item))
        push!(item_mat["Product code"], get_code(get_product(item)))
        push!(item_mat["Package code"], error("What is a package code?")) # looks like a little label to use for visualisation
        push!(item_mat["Number of items"], error("What does that mean?"))
        # An item with an item ident can exist in several copies. TODO
        push!(item_mat["Length"], get_dim(item).le)
        push!(item_mat["Width"], get_dim(item).wi)
        push!(item_mat["Height"], get_height(item))
        push!(item_mat["Weight"], get_weight(item))
        push!(item_mat["Nesting height"], get_nesting_height(item))
        push!(item_mat["Stackability code"], get_stackability_code(item))
        push!(item_mat["Forced orientation"], get_forced_orientation(item))
        push!(item_mat["Earliest arrival time,"], get_time_window(item)[1])
        push!(item_mat["Latest arrival time"], get_time_window(item)[2])
        push!(item_mat["Inventory cost"], get_inventory_cost(item))
        push!(item_mat["Max stackability"], get_max_stackability(get_product(item)))

    end

    # write truck infos
    # Supplier code, Supplier loading order, Supplier dock, Supplier dock loading order,
    # Plant code, Plant dock, Plant dock loading order, Product code, Arrival time,
    # Id truck, Length, Width, Height, Max weight, Stack with multiple docks,
    # Max density, Max weight on the bottom item in stacks, Cost, EMmm, EMmr, CM
    # CJfm, CJfc, CJfh, EM, EJhr, EJcr, EJeh
    

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