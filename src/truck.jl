include("dim.jl")
include("pos.jl")

struct Truck
    dim::Dim
    height::Real
    max_stack_density::Real
    max_stack_weight::Real
    supplier_orders::Dict{String, Integer}
    supplier_dock_orders::Dict{String, Dict{String, Integer}}
    plant_dock_orders::Dict{String, Integer}

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

get_dim(truck::Truck) = truck.dim
get_height(truck::Truck) = truck.height
get_max_stack_density(truck::Truck) = truck.max_stack_density
get_max_stack_weight(truck::Truck) = truck.max_stack_weight
get_supplier_orders(truck::Truck) = truck.supplier_orders
get_supplier_dock_orders(truck::Truck) = truck.supplier_dock_orders
get_supplier_dock_order(truck::Truck, supplier, supplier_dock) = truck.supplier_dock_orders[supplier][supplier_dock]

get_plant_dock_orders(truck::Truck) = truck.plant_dock_orders

get_supplier_order(truck::Truck, supplier) = get_supplier_orders(truck)[supplier]
get_plant_dock_order(truck::Truck, plant_dock) = get_plant_dock_orders(truck)[plant_dock]

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

get_suppliers(truck::Truck) = collect(keys(get_supplier_orders(truck)))

function set_supplier_orders!(truck, d)
    merge!(truck.supplier_orders, d)
end
function set_supplier_dock_orders!(truck, d)
    merge!(truck.supplier_dock_orders, d)

end
function set_plant_dock_orders!(truck, d)
    merge!(truck.plant_dock_orders, d)

end
