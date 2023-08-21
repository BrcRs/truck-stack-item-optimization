include("placement.jl")

struct OrderedStack <: AbstractStack
    stack::Stack
    supplier_order::Integer
    supplier_dock_order::Integer
    plant_dock_order::Integer
end

get_dim(s::OrderedStack) = get_dim(s.stack)
get_pos(s::OrderedStack) = get_pos(s.stack)
supplier_order(s::OrderedStack) = s.supplier_order
supplier_dock_order(s::OrderedStack) = s.supplier_dock_order
plant_dock_order(s::OrderedStack) = s.plant_dock_order


function order_instance(instance::Dict{T, Stack}) where T <: Integer

    # sort stacks by x value
    sorted_instance = sort(collect(instance), by=p -> get_pos(p[2]).x)

    res = Vector{Pair{Integer, OrderedStack}}()

    supplier = 1
    supplier_dock = 1
    plant_dock = 1
    for p in sorted_instance
        push!(res, p[1] => OrderedStack(p[2], supplier, supplier_dock, plant_dock))
        if rand() < 0.5
            plant_dock += 1
        elseif rand() < 0.5
            supplier_dock += 1
            plant_dock = 1
        elseif rand() < 0.5
            supplier += 1
            supplier_dock = 1
            plant_dock = 1
        end
    end

    return res

end