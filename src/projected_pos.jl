include("placement.jl")

is_projected(pos::AbstractPos) = typeof(pos) == ProjectedPos


function dummy_dim(pos::ProjectedPos; precision=3)
    projdim = missing
    if get_orientation(pos) == :Vertical
        projdim = Dim(10.0^-precision, max(10.0^-precision, get_origin(pos).y - get_pos(pos).y))
    elseif get_orientation(pos) == :Horizontal
        projdim = Dim(max(10.0^-precision, get_origin(pos).x - get_pos(pos).x), 10.0^-precision)
    end
    return projdim
end

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


function is_intersected(pos::ProjectedPos, s::AbstractStack; precision=3, verbose=false)

    
    projdim = dummy_dim(pos; precision=precision)
    
    if verbose
        println()
        println("Overlap Y: ", overlapY(get_pos(pos), projdim, get_pos(s), get_dim(s); precision=precision))
        println("Overlap X: ", overlapX(get_pos(pos), projdim, get_pos(s), get_dim(s); precision=precision))
    end
    return  is_intersected(get_pos(pos), projdim, get_pos(s), get_dim(s); precision=precision, verbose=verbose)
end

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
                    println(Pos(get_pos(c).x, get_pos(p).y))
                end
                push!(to_add, Pos(get_pos(c).x, get_pos(p).y))
            elseif get_orientation(c) == :Horizontal
                if verbose
                    println(Pos(get_pos(p).x, get_pos(c).y))
                end
                push!(to_add, Pos(get_pos(p).x, get_pos(c).y))
            else
                error("Unknown orientation: $(get_orientation(c))")
            end
        end
    end
    return
end
