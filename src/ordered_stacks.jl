include("placement.jl")

abstract type AbstractOrderedStack <: AbstractStack end


struct OrderedStack <: AbstractOrderedStack
    stack::Stack
    supplier_order::Integer
    supplier_dock_order::Integer
    plant_dock_order::Integer
end

function OrderedStack(pos::Pos, dim::Dim, supplier_order, supplier_dock_order, plant_dock_order) 
    OrderedStack(Stack(pos, dim), supplier_order, supplier_dock_order, plant_dock_order)
end

get_dim(s::OrderedStack) = get_dim(s.stack)
get_pos(s::OrderedStack) = get_pos(s.stack)
supplier_order(s::OrderedStack) = s.supplier_order
supplier_dock_order(s::OrderedStack) = s.supplier_dock_order
plant_dock_order(s::OrderedStack) = s.plant_dock_order

function can_be_placed(solution, o, s::OrderedStack, W, orientation::Symbol; precision=3, verbose=false)
    condition = missing
    # Make sure no other stack already placed and with greater x position has a lower loading order
    # take the stack with lower x among those with x > o.x
    greater_xs = [p[2] for p in solution if get_pos(p[2]).x > o.x]
    if verbose
        println("greater_xs")
        display(greater_xs)
    end
    if !isempty(greater_xs)
        sorted_greater_xs = sort(greater_xs, by= y -> get_pos(y).x)
        min_x = get_pos(sorted_greater_xs[begin]).x
        mingreater_xs = filter(stack -> eqtol(get_pos(stack).x, min_x, 3), sorted_greater_xs)
        condition = true
        for stack in mingreater_xs
            condition *= leq_order(s, stack)
        end
        if verbose
            println("mingreater_x")
            display(mingreater_x)
            println("condition: $condition")
        end
    else
        condition = true
    end
    return condition && can_be_placed(solution, o, s.stack, W, orientation; precision)

end

"""
Return true if s1 is before s2 or equal priority.
"""
function leq_order(s1, s2)

    if s1.supplier_order < s2.supplier_order
        return true
    elseif s1.supplier_order == s2.supplier_order && s1.supplier_dock_order < s2.supplier_dock_order
        return true
    elseif s1.supplier_order == s2.supplier_order && s1.supplier_dock_order == s2.supplier_dock_order && s1.plant_dock_order < s2.plant_dock_order
        return true
    elseif s1.supplier_order == s2.supplier_order && s1.supplier_dock_order == s2.supplier_dock_order && s1.plant_dock_order == s2.plant_dock_order
        return true
    else
        return false
    end
end

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