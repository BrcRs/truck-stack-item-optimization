include("placement.jl")

@auto_hash_equals struct IntersectionPos <: AbstractPos
    x
    y
end

function readable(p::IntersectionPos)
    return string("IntersectionPos(", round(p.x, digits=3),", ", round(p.y, digits=3),")")
end
get_pos(pos::IntersectionPos) = Pos(pos.x, pos.y)
# @auto_hash_equals struct Pos <: AbstractPos
#     x
#     y
# end

# function readable(p::Pos)
#     return string("Pos(", round(p.x, digits=3),", ", round(p.y, digits=3),")")
# end

"""
A ProjectedPos is a position whose origin is not against the ground nor the driver's 
cabin on the left. Thus, we imagine a line going from the origin to the left if 
Horizontal or to the ground if Vertical. The actual position (which is "projected") 
is the most distant position from the origin on the line without anything between it and the origin. 
"""
struct ProjectedPos <: AbstractPos
    p::Base.RefValue{Pos}
    origin::Pos
    orientation::Symbol
    function ProjectedPos(p, origin, orientation)
        if !(orientation in [:Horizontal, :Vertical])
            throw(ArgumentError("orientation field should be either :Horizontal or :Vertical.\nGot $orientation"))
        end
        if orientation == :Vertical && get_pos(p).x != get_pos(origin).x
            throw(ArgumentError("Vertical projected position and origin don't have same x.\n$p\n$origin"))
        elseif orientation == :Horizontal && get_pos(p).y != get_pos(origin).y
            throw(ArgumentError("Horizontal projected position and origin don't have same y.\n$p\n$origin"))
        end
        new(Ref(p), origin, orientation)
    end

end


Base.hash(a::ProjectedPos, h::UInt) = hash(a.p[], hash(a.origin, hash(a.orientiation, hash(:ProjectedPos, h))))
Base.:(==)(a::ProjectedPos, b::ProjectedPos) = isequal(a.p[], b.p[]) && isequal(a.origin, b.origin) && isequal(a.orientation, b.orientation)


# function Base.isequal(a::ProjectedPos, b::ProjectedPos)
#     return get_pos(a) == get_pos(b) && get_origin(a) == get_origin(b) && get_orientation(a) == get_orientation(b)
# end

# function Base.is(a::ProjectedPos, b::ProjectedPos)
#     error("Identity comparison between ProjectedPos not implemented")
# end

get_pos(pos::ProjectedPos) = pos.p[]
get_origin(pos::ProjectedPos) = pos.origin
get_orientation(pos::ProjectedPos) = pos.orientation

function set_pos!(propos::ProjectedPos, pos::Pos)
    # get_pos(propos) = pos
    propos.p[] = pos
end

is_projected(pos::AbstractPos) = typeof(pos) == ProjectedPos
is_intersectionPos(pos::AbstractPos) = typeof(pos) == IntersectionPos

is_secure_pos(pos::AbstractPos) = !is_projected(pos) && !is_intersectionPos(pos)

function dummy_dim(pos::ProjectedPos; precision=3)
    projdim = missing
    if get_orientation(pos) == :Vertical
        projdim = Dim(10.0^-precision, max(10.0^-precision, get_origin(pos).y - get_pos(pos).y))
    elseif get_orientation(pos) == :Horizontal
        projdim = Dim(max(10.0^-precision, get_origin(pos).x - get_pos(pos).x), 10.0^-precision)
    end
    return projdim
end


"""Given a position, find the most to the left available position without 
overlapping a placed stack."""
function totheleft(pos::Pos, solution; precision=3)
    # Find all stacks overlapping on y axis and with x < pos.x
    boxesleft = findboxesleft(pos, Dim(10.0^-precision, 10.0^-precision), solution; precision=precision)

    # Find the stack which extends the most to the right
    rightsides = [get_pos(solution[k]).x + get_dim(solution[k]).le for k in boxesleft]
    leftbound = isempty(rightsides) ? 0 : max(rightsides...)    

    # return the Pos with x position as the right side of the stack
    return eqtol(leftbound, get_pos(pos).x, precision) ? 
                    Pos(get_pos(pos).x, get_pos(pos).y) : 
       ProjectedPos(Pos(leftbound,      get_pos(pos).y), get_pos(pos), :Horizontal)
end


"""Given a position, find the most to the bottom available position without 
overlapping a placed stack."""
function tothebottom(pos::Pos, solution; precision=3)
    # Find all stacks overlapping on x axis and with y < pos.y
    boxesbot = findboxesbelow(pos, Dim(10.0^-precision, 10.0^-precision), solution; precision)

    # Find the stack which extends the most to the top
    topsides = [get_pos(solution[k]).y + get_dim(solution[k]).wi for k in boxesbot]
    botbound = isempty(topsides) ? 0 : max(topsides...)    

    # return the Pos with y position as the top side of the stack
    return eqtol(botbound, pos.y, precision) ? 
                    Pos(pos.x, pos.y) :
       ProjectedPos(Pos(pos.x, botbound), pos, :Vertical)
end


"""
    is_intersected(pos::ProjectedPos, s::AbstractStack; precision=3, verbose=false)

Return if ProjectedPos pos is intersected by s, or in other terms, if s is 
between pos's origin and pos's projected pos.
"""
function is_intersected(pos::ProjectedPos, s::AbstractStack; precision=3, verbose=false)

    
    projdim = dummy_dim(pos; precision=precision)
    
    if verbose
        println()
        println("Overlap Y: ", overlapY(get_pos(pos), projdim, get_pos(s), get_dim(s); precision=precision))
        println("Overlap X: ", overlapX(get_pos(pos), projdim, get_pos(s), get_dim(s); precision=precision))
    end
    return  is_intersected(get_pos(pos), projdim, get_pos(s), get_dim(s); precision=precision, verbose=verbose)
end

"""
    is_intersected(p1::ProjectedPos, p2::ProjectedPos; precision=3, verbose=false)

Return if ProjectedPos p1 is intersected by p2::ProjectedPos, or in other terms, if p2's line is 
between p1's origin and p1's projected pos.
"""
function is_intersected(p1::ProjectedPos, p2::ProjectedPos; precision=3, verbose=false)

    projdim1 = dummy_dim(p1; precision=precision)

    projdim2 = dummy_dim(p2; precision=precision)
    
    # if verbose
    #     println()
    #     println("Overlap Y: ", overlapY(get_pos(pos), projdim, get_pos(s), get_dim(s); precision=precision))
    #     println("Overlap X: ", overlapX(get_pos(pos), projdim, get_pos(s), get_dim(s); precision=precision))
    # end
    return  is_intersected(get_pos(p1), projdim1, get_pos(p2), projdim2; precision=precision, verbose=verbose)
end

"""

Update the projected position of o::ProjectedPos in function of intersecting s::AbstractStack.
If the stack covers the origin, o is removed. Else, we find the most distant 
position on the line of o without breaking the line (this position is given by s, 
it is either the y value of its top or the x value of its right side).
"""
function upd!(corners, o::ProjectedPos, s::AbstractStack; precision=3, verbose=false)
    if get_orientation(o) == :Vertical
        # if the top of the stack is higher than the origin, the whole projected
        # corner is covered, thus we must delete it

        if verbose
            println("Stack $s")
            println("greatertol(get_pos(s).y + get_dim(s).wi, get_origin(o).y, precision)")
            println("greatertol($(get_pos(s).y) + $(get_dim(s).wi), $(get_origin(o).y), $precision) = ", greatertol(get_pos(s).y + get_dim(s).wi, get_origin(o).y, precision))
            println()
        end

        if greatertol(get_pos(s).y + get_dim(s).wi, get_origin(o).y, precision)
            filter!(x -> x != o, corners)
        else
            # else, give it the value of the height of the top of the stack
            set_pos!(o, Pos(get_pos(o).x, get_pos(s).y + get_dim(s).wi))
        end

    elseif get_orientation(o) == :Horizontal

        if greatertol(get_pos(s).x + get_dim(s).le, get_origin(o).x, precision)
            filter!(x -> x != o, corners)
        else
            set_pos!(o, Pos(get_pos(s).x + get_dim(s).le, get_pos(o).y))
        end

    else
        throw(ArgumentError("Unknown orientation \"$orientation\""))
    end

    return
end

"""
    upd_intersection!(to_add, c::ProjectedPos, allprojected::Vector{ProjectedPos}; verbose=false)

c is a new ProjectedPos. Iterate over all known projectedpos to determine is c 
intersect any. If yes, then create a new position at the intersection (added to `to_add`).
"""
function upd_intersection!(to_add, c::ProjectedPos, allprojected::Vector{ProjectedPos}; verbose=false)
    # for each projected pos
    for p in allprojected
        # check if orthogonal
        # check if there is intersection with c
        if verbose
            println(p)
            println("get_orientation($c) != get_orientation($p) ", get_orientation(c) != get_orientation(p))
            println("s_intersected($c, $p) ", is_intersected(c, p))
        end
        if get_orientation(c) != get_orientation(p) && is_intersected(c, p)
        # if yes, add a new corner to to_add at the intersection
            if get_orientation(c) == :Vertical
                if verbose
                    println(IntersectionPos(get_pos(c).x, get_pos(p).y))
                end
                push!(to_add, IntersectionPos(get_pos(c).x, get_pos(p).y))
            elseif get_orientation(c) == :Horizontal
                if verbose
                    println(IntersectionPos(get_pos(p).x, get_pos(c).y))
                end
                push!(to_add, IntersectionPos(get_pos(p).x, get_pos(c).y))
            else
                error("Unknown orientation: $(get_orientation(c))")
            end
        end
    end
    return
end
