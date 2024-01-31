using Random
using AutoHashEquals
# 1 => length
# 2 => width

include("truck.jl")
include("dim.jl")
include("pos.jl")

readable(a::Nothing) = "nothing"

# function Dim(le, wi)
#     if le == 0 || wi == 0
#         throw(ArgumentError("Dimension can't be of size zero"))
#     end
#     return Dim(le, wi)
# end


abstract type AbstractStack end
struct Stack <: AbstractStack
    pos::Pos
    dim::Dim
end

get_pos(pos::Pos) = pos


function newStack(s::Stack, p::Pos, d::Dim)
    return Stack(p, d)
end



function is_intersected(pos1::Pos, dim1::Dim, pos2::Pos, dim2::Dim; precision=3, verbose=false)
    
    return  overlapY(pos1, dim1, pos2, dim2; precision=precision) &&
            overlapX(pos1, dim1, pos2, dim2; precision=precision)
end


get_dim(s::Stack) = s.dim
get_pos(s::Stack) = s.pos



"""
    findboxesabove(pos::Pos, r::Dict{T, S}; precision=3) where {T <: Integer, S <: AbstractStack}

Return elements of `r` that are directly above and overlapping on x axis with `pos`.
By convention, the right and top borders of a stack aren't counted in the stack.
"""
function findboxesabove(pos::Pos, r::Dict{T, S}; precision=3) where {T <: Integer, S <: AbstractStack}
    # return filter(b -> get_pos(r[b]).y >= pos.y && get_pos(r[b]).x < pos.x + dim.le && get_pos(r[b]).x + get_dim(r[b]).le > pos.x, keys(r))
    # return filter(b -> geqtol(get_pos(r[b]).y, pos.y, precision) && lessertol(get_pos(r[b]).x, pos.x + dim.le, precision) && greatertol(get_pos(r[b]).x + get_dim(r[b]).le, pos.x, precision), keys(r))
    return filter(b -> geqtol(get_pos(r[b]).y, pos.y, precision) && leqtol(get_pos(r[b]).x, pos.x, precision) && greatertol(get_pos(r[b]).x + get_dim(r[b]).le, pos.x, precision), keys(r))
end

"""
    findboxesbelow(pos::Pos, dim::Dim, r::Dict{T, S}; precision=3) where {T <: Integer, S <: AbstractStack}

Return elements of `r` that are directly below and overlapping on x axis with a 
stack of position `pos` and dimensions `dim`.
By convention, the right and top borders of a stack aren't counted in the stack.
"""
function findboxesbelow(pos::Pos, dim::Dim, r::Dict{T, S}; precision=3) where {T <: Integer, S <: AbstractStack}
    # return filter(b -> get_pos(r[b]).y >= pos.y && get_pos(r[b]).x < pos.x + dim.le && get_pos(r[b]).x + get_dim(r[b]).le > pos.x, keys(r))
    # return filter(b -> geqtol(get_pos(r[b]).y, pos.y, precision) && lessertol(get_pos(r[b]).x, pos.x + dim.le, precision) && greatertol(get_pos(r[b]).x + get_dim(r[b]).le, pos.x, precision), keys(r))
    return filter(b ->  leqtol(get_pos(r[b]).y, pos.y, precision) && 
                        leqtol(get_pos(r[b]).x, pos.x + dim.le, precision) && 
                        greatertol(get_pos(r[b]).x + get_dim(r[b]).le, pos.x, precision), keys(r))
end

"""
    findboxesright(pos::Pos, r::Dict{T, S}; precision::Integer=3) where {T <: Integer, S <: AbstractStack}

Return elements of `r` that are directly to the right and overlapping on y axis with `pos`.
By convention, the right and top borders of a stack aren't counted in the stack.
"""
function findboxesright(pos::Pos, r::Dict{T, S}; precision::Integer=3) where {T <: Integer, S <: AbstractStack}
    # return filter(b -> geqtol(get_pos(r[b]).x, pos.x, precision) && lessertol(get_pos(r[b]).y, pos.y + dim.wi, precision) && greatertol(get_pos(r[b]).y + get_dim(r[b]).wi, pos.y, precision), keys(r))
    return filter(
        b -> geqtol(get_pos(r[b]).x, pos.x, precision) && 
        leqtol(get_pos(r[b]).y, pos.y, precision) && 
        greatertol(get_pos(r[b]).y + get_dim(r[b]).wi, pos.y, precision), 
        keys(r))
end

"""
    findboxesright(pos::Pos, dim::Dim, r::Dict{T, S}; precision=3) where {T <: Integer, S <: AbstractStack}

Return elements of `r` that are directly to the right and overlapping on y axis 
with a stack of position `pos` and dimensions `dim`.
By convention, the right and top borders of a stack aren't counted in the stack.
"""
function findboxesright(pos::Pos, dim::Dim, r::Dict{T, S}; precision=3) where {T <: Integer, S <: AbstractStack}
    return filter(
        b -> geqtol(get_pos(r[b]).x, pos.x, precision) && 
        lessertol(get_pos(r[b]).y, pos.y + dim.wi, precision) && 
        greatertol(get_pos(r[b]).y + get_dim(r[b]).wi, pos.y, precision), 
        keys(r))
    # return filter(
    #     b -> geqtol(get_pos(r[b]).x, pos.x, precision) && 
    #     leqtol(get_pos(r[b]).y, pos.y, precision) && 
    #     greatertol(get_pos(r[b]).y + get_dim(r[b]).wi, pos.y, precision), 
    #     keys(r))
end

"""
    findboxesleft(pos::Pos, dim::Dim, r::Dict{T, S}; precision=3) where {T <: Integer, S <: AbstractStack}

Return elements of `r` that are directly to the left and overlapping on y axis 
with a stack of position `pos` and dimensions `dim`.
By convention, the right and top borders of a stack aren't counted in the stack.
"""
function findboxesleft(pos::Pos, dim::Dim, r::Dict{T, S}; precision=3, verbose=false) where {T <: Integer, S <: AbstractStack}
    # return filter(
    #     b -> leqtol(get_pos(r[b]).x, pos.x, precision) && 
    #     lessertol(get_pos(r[b]).y, pos.y + dim.wi, precision) && 
    #     greatertol(get_pos(r[b]).y + get_dim(r[b]).wi, pos.y, precision),
    #     keys(r))
    if verbose
        for b in keys(r)
            println(b)
            println("leqtol(get_pos(r[b]).x, pos.x, precision)\n", leqtol(get_pos(r[b]).x, pos.x, precision))
            println("overlapY(pos, dim, get_pos(r[b]), get_dim(r[b]); precision=precision)\n", overlapY(pos, dim, get_pos(r[b]), get_dim(r[b]); precision=precision))
        end
    end

    return filter(
        b -> leqtol(get_pos(r[b]).x, pos.x, precision) && 
        overlapY(pos, dim, get_pos(r[b]), get_dim(r[b]); precision=precision),
        keys(r))
end

"""Return wether two boxes overlap on X"""
function overlapX(apos, adim, bpos, bdim; precision=3)
    return greatertol(apos.x + adim.le, bpos.x, precision) && lessertol(apos.x, bpos.x + bdim.le, precision)
end


"""Return wether two boxes overlap on Y"""
function overlapY(apos, adim, bpos, bdim; precision=3)
    return greatertol(apos.y + adim.wi, bpos.y, precision) && lessertol(apos.y, bpos.y + bdim.wi, precision)
end

"""<= but with tolerance"""
leqtol(a, b, decimals=3) = round(a, digits=decimals) <= round(b, digits=decimals)

"""> but with tolerance"""
greatertol(a, b, decimals=3) = !eqtol(a, b, decimals) && geqtol(a, b, decimals)

"""< but with tolerance"""
lessertol(a, b, decimals=3) = !eqtol(a, b, decimals) && leqtol(a, b, decimals)

""">= but with tolerance"""
geqtol(a, b, decimals=3) = round(a, digits=decimals) >= round(b, digits=decimals)

"""== but with tolerance"""
eqtol(a, b, decimals=3) = round(a, digits=decimals) == round(b, digits=decimals)

"""
    collision(pos::Pos, dim::Dim, r::Dict{T, S}; precision=3, verbose=false) where {T <: Integer, S <: AbstractStack}

Return true if considered box overlaps with existing one in solution `r`.
"""
function collision(pos::Pos, dim::Dim, r::Dict{T, S}; precision=3, verbose=false) where {T <: Integer, S <: AbstractStack}

    # Find boxes with corresponding x
    # tocheck = filter(b -> get_pos(r[b]).x < pos.x + dim.le && get_pos(r[b]).x + get_dim(r[b]).le > pos.x, keys(r))
    # tocheck = findboxesabove(pos, dim, r, precision)
    # check each one
    for k in keys(r)

        if overlapX(pos, dim, get_pos(r[k]), get_dim(r[k]); precision) && overlapY(pos, dim, get_pos(r[k]), get_dim(r[k]); precision)
            if verbose
                println(k)
                println("overlapX($pos, $dim, $(get_pos(r[k])), $(get_dim(r[k])); $precision) = ", overlapX(pos, dim, get_pos(r[k]), get_dim(r[k]); precision))
                println("overlapY($pos, $dim, $(get_pos(r[k])), $(get_dim(r[k])); $precision) = ", overlapY(pos, dim, get_pos(r[k]), get_dim(r[k]); precision))
            end
            return true
        end
    end

    return false
end

function outofbound(pos, dim, W; precision=3)
    return greatertol(pos.y + dim.wi, W, precision)
end

"""
    coveredcorners(corners::Vector{<:AbstractPos}, o::Pos, le, wi; precision=3, verbose=false)

Remove corners covered by provided stack at position o.
"""
function coveredcorners(corners::Vector{<:AbstractPos}, o::Pos, le, wi; precision=3, verbose=false)
    torem = AbstractPos[]
    for o2 in corners
        if leqtol(o.x, get_pos(o2).x, precision) && lessertol(get_pos(o2).x, o.x + le, precision) && leqtol(o.y, get_pos(o2).y, precision) && lessertol(get_pos(o2).y, o.y + wi, precision)
            push!(torem, o2)
            if verbose
                println("-==-")
                println("leqtol($(o.x), $(get_pos(o2).x), $precision) = ", leqtol(o.x, get_pos(o2).x, precision) )
                println("lessertol($(get_pos(o2).x), $(o.x) + $(le), $precision) = ", lessertol(get_pos(o2).x, o.x + le, precision) )
                println("leqtol($(o.y), $(get_pos(o2).y), $precision) = ", leqtol(o.y, get_pos(o2).y, precision) )
                println("lessertol($(get_pos(o2).y), $(o.y) + $(wi), $precision) = ", lessertol(get_pos(o2).y, o.y + wi, precision))
                println("-==-")
            end
        end
    
    end
    return torem
end


"""
    is_secure(stack, solution; precision=3, verbose=false)

Return `true` if the stack is either against the truck's driver cabin or directly against (on its left) to another stack.  
"""
function is_secure(stack, solution; precision=3, verbose=false)
    if eqtol(get_pos(stack).x, 0, precision)
        if verbose
            println("The stack is against the cabin.")
        end
        return true
    end
    boxesleft = findboxesleft(get_pos(stack), get_dim(stack), solution; precision=precision, verbose=verbose)
    if verbose
        println("is_secure examines all boxes to the left")
    end
    for k in boxesleft
        b = solution[k]
        if verbose
            println(k)
            println("eqtol(get_pos(b).x + get_dim(b).le, get_pos(stack).x, precision)\n", eqtol(get_pos(b).x + get_dim(b).le, get_pos(stack).x, precision))
        end
        if eqtol(get_pos(b).x + get_dim(b).le, get_pos(stack).x, precision)
            if verbose
                println("The stack is against stack $k.")
            end
            return true
        end
    end
    if verbose
        println("There isn't any stack the stack is against, and is not against the cabin.")
    end
    return false
end

function can_be_placed(solution, o::Pos, s::Stack, W, orientation::Symbol; precision=3, verbose=false)
    return can_be_placed(solution, o, get_dim(s), W, orientation; precision=precision, verbose=verbose)
end

"""
    can_be_placed(solution, o::Pos, dim::Dim, truck::Truck, orientation::Symbol; precision=3, verbose=false)

Return `true` if a stack of dimension dim can be placed at position o without colliding with another stack or getting out of bounds.
"""
function can_be_placed(solution, o::Pos, dim::Dim, truck::Truck, orientation::Symbol; precision=3, verbose=false)
    W = get_dim(truck).wi
    res = nothing

    if orientation == :Perpendicular
        res = leqtol(o.y + dim.le, W, precision) && !collision(Pos(o.x, o.y), Dim(dim.wi, dim.le), solution; precision)
    elseif orientation == :Parallel
        res = leqtol(o.y + dim.wi, W, precision) && !collision(Pos(o.x, o.y), Dim(dim.le, dim.wi), solution; precision)
    else
        throw(ArgumentError("orientation argument should either be :Perpendicular or :Parallel.\nUnknown flag: $orientation"))
    end 
    return res
    
end


"""Return wether position o is illegal with already placed boxes."""
function illegalpos(o, r; precision=3)
    return collision(o, Dim(1/(10^precision), 1/(10^precision)), r)
end


