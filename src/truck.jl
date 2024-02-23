include("dim.jl")
include("pos.jl")
include("product.jl")

struct Truck
    id::String
    dim::Dim
    height::Real
    max_stack_density::Real
    max_stack_weights::Dict{String, Real} # string is product code
    TMm::Real # Max authorized loading weight of truck t
    supplier_orders::Dict{Integer, Integer}
    supplier_dock_orders::Dict{Integer, Dict{String, Integer}}
    plant_dock_orders::Dict{String, Integer}

    arrival_time::Any
    cost::Integer

    CM::Real # weight of the tractor
    CJ_fm::Real #  distance between the front and middle axles of the tractor
    CJ_fc::Real # distance between the front axle and the center of gravity of the tractor
    CJ_fh::Real # distance between the front axle and the harness of the tractor
    EM::Real #  weight of the empty trailer
    EJ_hr::Real # distance between the harness and the rear axle of the trailer
    EJ_cr::Real #  distance between the center of gravity of the trailer and the rear axle
    EJ_eh::Real #  distance between the start of the trailer and the harness
    EM_mr::Real # max weight on the rear axle of the trailer
    EM_mm::Real # max weight on the middle axle of the trailer
end

function Truck(
    dim::Dim, height, max_stack_density, max_stack_weights, TMm, supplier_orders, 
    supplier_dock_orders, plant_dock_orders, arrival_time, cost, CM, CJ_fm, CJ_fc, CJ_fh, EM,
    EJ_hr, EJ_cr, EJ_eh, EM_mr, EM_mm
) 
    return Truck("", dim, height, max_stack_density, max_stack_weights, TMm, supplier_orders, 
    supplier_dock_orders, plant_dock_orders, arrival_time, cost, CM, CJ_fm, CJ_fc, CJ_fh, EM,
    EJ_hr, EJ_cr, EJ_eh, EM_mr, EM_mm)
end

function add_max_stack_weights!(truck::Truck, product::Product, max_stack_weight)
    truck.max_stack_weights[get_code(product)] = max_stack_weight
end

get_id(truck::Truck) = truck.id
get_dim(truck::Truck) = truck.dim
get_height(truck::Truck) = truck.height
get_max_stack_density(truck::Truck) = truck.max_stack_density
get_max_stack_weights(truck::Truck) = truck.max_stack_weights
get_supplier_orders(truck::Truck) = truck.supplier_orders
get_supplier_dock_orders(truck::Truck) = truck.supplier_dock_orders
get_supplier_dock_order(truck::Truck, supplier, supplier_dock) = truck.supplier_dock_orders[supplier][supplier_dock]

get_plant_dock_orders(truck::Truck) = truck.plant_dock_orders

get_supplier_order(truck::Truck, supplier) = get_supplier_orders(truck)[supplier]
get_plant_dock_order(truck::Truck, plant_dock) = get_plant_dock_orders(truck)[plant_dock]
get_arrival_time(truck::Truck) = truck.arrival_time
get_cost(truck::Truck) = truck.cost

get_CM(truck::Truck) = truck.CM
get_CJ_fm(truck::Truck) = truck.CJ_fm
get_CJ_fc(truck::Truck) = truck.CJ_fc
get_CJ_fh(truck::Truck) = truck.CJ_fh
get_EM(truck::Truck) = truck.EM
get_EJ_hr(truck::Truck) = truck.EJ_hr
get_EJ_cr(truck::Truck) = truck.EJ_cr
get_EJ_eh(truck::Truck) = truck.EJ_eh
get_EM_mr(truck::Truck) = truck.EM_mr
get_EM_mm(truck::Truck) = truck.EM_mm

get_TMm(truck::Truck) = truck.TMm

get_suppliers(truck::Truck) = collect(keys(get_supplier_orders(truck)))

get_volume(truck::Truck) = get_area(truck) * get_height(truck)
get_area(truck::Truck) = get_dim(truck).le * get_dim(truck).wi

function set_supplier_orders!(truck, d)
    merge!(truck.supplier_orders, d)
end
function set_supplier_dock_orders!(truck, d)
    merge!(truck.supplier_dock_orders, d)

end
function set_plant_dock_orders!(truck, d)
    merge!(truck.plant_dock_orders, d)

end


function set_height(truck::Truck, h::Real)
    return Truck(
        truck.id,
        truck.dim,
        h,
        truck.max_stack_density,
        truck.max_stack_weights,
        truck.TMm,
        truck.supplier_orders,
        truck.supplier_dock_orders,
        truck.plant_dock_orders,
        truck.arrival_time,
        truck.cost,

        truck.CM,
        truck.CJ_fm,
        truck.CJ_fc,
        truck.CJ_fh,
        truck.EM,
        truck.EJ_hr,
        truck.EJ_cr,
        truck.EJ_eh,
        truck.EM_mr,
        truck.EM_mm,
    )
end
function set_id(truck::Truck, id::String)
    return Truck(
        id,
        truck.dim,
        truck.height,
        truck.max_stack_density,
        truck.max_stack_weights,
        truck.TMm,
        truck.supplier_orders,
        truck.supplier_dock_orders,
        truck.plant_dock_orders,
        truck.arrival_time,
        truck.cost,

        truck.CM,
        truck.CJ_fm,
        truck.CJ_fc,
        truck.CJ_fh,
        truck.EM,
        truck.EJ_hr,
        truck.EJ_cr,
        truck.EJ_eh,
        truck.EM_mr,
        truck.EM_mm,
    )
end