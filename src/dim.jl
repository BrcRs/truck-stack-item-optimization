using AutoHashEquals


@auto_hash_equals struct Dim
    le # le:length corresponds to x axis of a truck
    wi # wi:width corresponds to y axis
    Dim(a, b) = a == 0 || b == 0 ? throw(ArgumentError("Dim($a, $b): Dimension can't be of size zero")) : new(a, b)
end

function readable(d::Dim)
    return string("Dim(", round(d.le, digits=3),", ", round(d.wi, digits=3),")")
end