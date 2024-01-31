
abstract type AbstractPos end

@auto_hash_equals struct Pos <: AbstractPos
    x
    y
end

function readable(p::Pos)
    return string("Pos(", round(p.x, digits=3),", ", round(p.y, digits=3),")")
end