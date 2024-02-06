using CSV
using Tables
using OrderedCollections

include("item.jl")



function return_label(i)
    # generator-like
    # A, B, C... etc then BA, BB, BC, etc...
    letters = 'A':'Z'
    return return_label(letters, i, "")
end

function return_label(characters, i, res)
    if i <= length(characters)
        return characters[i] * res
    else
        return return_label(
            characters, 
            div(i, length(characters)),
            characters[i % length(characters) + 1] * res
        )
    end

end

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
        stack_label = return_label(i)
        push!(stack_mat["Stack code"], stack_label) # DONE stack code is a simple label to easily identify the stack
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
            push!(item_mat["Item code"], string(stack_label, get_copy_number(item)))
            # error("Item code is stack code + item copy number")
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



    # DONE convert to table
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

    # DONE then write to directory
    # output_items.csv, output_stacks.csv and output_trucks.csv
    CSV.write(joinpath(directory, "output_trucks.csv"), truck_tab, append=append, writeheader=!append, delim=";", decimal=',')
    CSV.write(joinpath(directory, "output_stacks.csv"), stack_tab, append=append, writeheader=!append, delim=";", decimal=',')
    CSV.write(joinpath(directory, "output_items.csv"), item_tab, append=append, writeheader=!append, delim=";", decimal=',')

end

function write_input(truck, items, directory; append=true)

    truck_mat = OrderedDict{String, Vector{Any}}()
    item_mat = OrderedDict{String, Vector{Any}}()


    truck_cols = [
        "Supplier code", "Supplier loading order", "Supplier dock", "Supplier dock loading order",
        "Plant code", "Plant dock", "Plant dock loading order", "Product code", "Arrival time",
        "Id truck", "Length", "Width", "Height", "Max weight", "Stack with multiple docks",
        "Max density", "Max weight on the bottom item in stacks", "Cost", "EMmm", "EMmr", "CM",
        "CJfm", "CJfc", "CJfh", "EM", "EJhr", "EJcr", "EJeh"
    ]
    for col in truck_cols
        truck_mat[col] = Vector{Any}()
    end

    item_cols = [
        "Item ident", "Supplier code", "Supplier dock", "Plant code", "Plant dock", 
        "Product code", "Package code", "Number of items", "Length", "Width", "Height", "Weight",
        "Nesting height", "Stackability code", "Forced orientation", "Earliest arrival time",
        "Latest arrival time", "Inventory cost", "Max stackability"
    ]
    for col in item_cols
        item_mat[col] = Vector{Any}()
    end

    # write item infos
    # Item ident, Supplier code, Supplier dock, Plant code, Plant dock, 
    # Product code, Package code, Number of items, Length, Width, Height, Weight,
    # Nesting height, Stackability code, Forced orientation, Earliest arrival time,
    # Latest arrival time, Inventory cost, Max stackability

    count_packages = Dict{String, Integer}()

    for item in items
        count_packages[get_package_code(item)] = 
            haskey(count_packages, get_package_code(item)) ? 
            count_packages[get_package_code(item)] + 1 : 1
    end

    products = Set{Product}()

    for item in items
        push!(products, get_product(item))

        push!(item_mat["Item ident"], get_id(item))
        push!(item_mat["Supplier code"], get_supplier(item))
        push!(item_mat["Supplier dock"], get_supplier_dock(item))
        push!(item_mat["Plant code"], get_plant(item))
        push!(item_mat["Plant dock"], get_plant_dock(item))
        push!(item_mat["Product code"], get_code(get_product(item)))
        push!(item_mat["Package code"], get_package_code(item)) # Copies of the same item share package code
        push!(item_mat["Number of items"], count_packages[get_package_code(item)])
        # An item with an item ident can exist in several copies. DONE
        push!(item_mat["Length"], get_dim(item).le)
        push!(item_mat["Width"], get_dim(item).wi)
        push!(item_mat["Height"], get_height(item))
        push!(item_mat["Weight"], get_weight(item))
        push!(item_mat["Nesting height"], get_nesting_height(item))
        push!(item_mat["Stackability code"], get_stackability_code(item))
        push!(item_mat["Forced orientation"], get_forced_orientation(item))
        push!(item_mat["Earliest arrival time"], get_time_window(item)[1])
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
    

    for supplier in keys(get_supplier_orders(truck))
        for supplier_dock in keys(get_supplier_dock_orders(truck)[supplier])
            for plant_dock in keys(get_plant_dock_orders(truck))
                for product in products
                    push!(truck_mat["Supplier code"], supplier)
                    push!(truck_mat["Supplier loading order"], get_supplier_orders(truck)[supplier])
                    push!(truck_mat["Supplier dock"], supplier_dock)
                    push!(truck_mat["Supplier dock loading order"], get_supplier_dock_order(truck, supplier, supplier_dock))
                    push!(truck_mat["Plant code"], 0) # TODO add plant code attribute to Truck
                    push!(truck_mat["Plant dock"], plant_dock)
                    push!(truck_mat["Plant dock loading order"], get_plant_dock_orders(truck)[plant_dock])
                    push!(truck_mat["Product code"], get_code(product))
                    push!(truck_mat["Arrival time"], 0)
                    push!(truck_mat["Id truck"], get_id(truck))
                    push!(truck_mat["Length"], get_dim(truck).le)
                    push!(truck_mat["Width"], get_dim(truck).wi)
                    push!(truck_mat["Height"], get_height(truck))
                    push!(truck_mat["Max weight"], get_TMm(truck))
                    # error("Continue implementing this function")
                    push!(truck_mat["Stack with multiple docks"], 0) # TODO implement this attribute in OrderedStack
                    push!(truck_mat["Max density"], get_max_stack_density(truck))
                    push!(truck_mat["Max weight on the bottom item in stacks"], get_max_stack_weights(truck)[get_code(product)])
                    push!(truck_mat["Cost"], get_cost(truck)) # DONE add cost to truck attributes 
                    push!(truck_mat["EMmm"], get_EM_mm(truck))
                    push!(truck_mat["EMmr"], get_EM_mr(truck))
                    push!(truck_mat["CM"], get_CM(truck))
                    push!(truck_mat["CJfm"], get_CJ_fm(truck))
                    push!(truck_mat["CJfc"], get_CJ_fc(truck))
                    push!(truck_mat["CJfh"], get_CJ_fh(truck))
                    push!(truck_mat["EM"], get_EM(truck))
                    push!(truck_mat["EJhr"], get_EJ_hr(truck))
                    push!(truck_mat["EJcr"], get_EJ_cr(truck))
                    push!(truck_mat["EJeh"], get_EJ_eh(truck))
                end
            end
        end
    end

    # DONE convert to table
    truck_tab = Tables.dictrowtable(truck_mat)
    item_tab = Tables.dictrowtable(item_mat)
    
    # display(truck_mat)
    # readline()
    # display(stack_mat)
    # readline()
    # display(item_mat)
    # readline()

    touch(joinpath(directory, "input_trucks.csv"))
    touch(joinpath(directory, "input_items.csv"))

    # DONE then write to directory
    # output_items.csv, output_stacks.csv and output_trucks.csv
    CSV.write(joinpath(directory, "input_trucks.csv"), truck_tab, append=append, writeheader=!append, delim=";", decimal=',')
    CSV.write(joinpath(directory, "input_items.csv"), item_tab, append=append, writeheader=!append, delim=";", decimal=',')
    

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