
@auto_hash_equals struct Product
    code::String
    max_stackability::Integer # max number of items an item of product can support
    max_weight::Real # max weight an item of product can support 
end

function Product(max_stackability, max_weight)
    return Product("", max_stackability, max_weight)
end

get_max_stackability(p::Product) = p.max_stackability
get_max_weight(p::Product) = p.max_weight
get_code(p::Product) = p.code
