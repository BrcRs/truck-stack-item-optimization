using Test
using Random

include("../src/item.jl")

@testset "rand_items" begin
    items = rand_items(
        100, 
        1, 
        10, 
        10, 
        100, 
        10, 
        100, 
        10, 
        randstring(8); 
        min_dim=0.001
        )

    display(items)

end


# @testset "make_stacks" begin
#     make_stacks(
#         items::Vector{Items}, 
#         plant_dock_orders, 
#         supplier_orders, 
#         supplier_dock_orders, 
#         max_height
#         )

# end