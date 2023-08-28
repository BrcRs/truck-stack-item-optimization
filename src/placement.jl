using Random
using AutoHashEquals
# 1 => length
# 2 => width

@auto_hash_equals struct Dim
    le # le:length corresponds to x axis of a truck
    wi # wi:width corresponds to y axis
    Dim(a, b) = a == 0 || b == 0 ? throw(ArgumentError("Dim($a, $b): Dimension can't be of size zero")) : new(a, b)
end

# function Dim(le, wi)
#     if le == 0 || wi == 0
#         throw(ArgumentError("Dimension can't be of size zero"))
#     end
#     return Dim(le, wi)
# end

abstract type AbstractPos end

@auto_hash_equals struct Pos <: AbstractPos
    x
    y
end

@auto_hash_equals struct ProjectedPos <: AbstractPos
    p::Pos
    origin::Pos
    orientation::Symbol
    function ProjectedPos(p, origin, orientation)
        if !(orientation in [:Horizontal, :Vertical])
            throw(ArgumentError("orientation field should be either :Horizontal or :Vertical.\nGot $orientation"))
        end
        new(p, origin, orientation)
    end

end

get_pos(pos::Pos) = pos
get_pos(pos::ProjectedPos) = pos.p
get_origin(pos::ProjectedPos) = pos.origin
get_orientation(pos::ProjectedPos) = pos.orientation

is_projected(pos::AbstractPos) = typeof(pos) == ProjectedPos

function is_intersected(pos::ProjectedPos, s::AbstractStack; precision=3)
    return  overlapY(get_pos(pos), Dim(10^-precision, 10^-precision), get_pos(s), get_dim(s); precision=precision) &&
            overlapX(get_pos(pos), Dim(10^-precision, 10^-precision), get_pos(s), get_dim(s); precision=precision) 
end

abstract type AbstractStack end
struct Stack <: AbstractStack
    pos::Pos
    dim::Dim
end

get_dim(s::Stack) = s.dim
get_pos(s::Stack) = s.pos

function upd!(corners, o::ProjectedPos, s::Stack)
    error("Not implemented yet")
end

# for o in corners
#     if greatertol(o.x, get_pos(solution[i]).x + get_dim(solution[i]).le)
#         break
#     end
#     if is_projected(o) && is_intersected(o, solution[i])
#         push!(to_upd, o)
#     end
# end
# for o in to_upd
#     upd!(corners, o, solution[i])
# end



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
    findboxesleft(pos::Pos, dim::Dim, r::Dict{T, S}, precision=3) where {T <: Integer, S <: AbstractStack}

Return elements of `r` that are directly to the left and overlapping on y axis 
with a stack of position `pos` and dimensions `dim`.
By convention, the right and top borders of a stack aren't counted in the stack.
"""
function findboxesleft(pos::Pos, dim::Dim, r::Dict{T, S}, precision=3) where {T <: Integer, S <: AbstractStack}
    return filter(
        b -> leqtol(get_pos(r[b]).x, pos.x, precision) && 
        lessertol(get_pos(r[b]).y, pos.y + dim.wi, precision) && 
        greatertol(get_pos(r[b]).y + get_dim(r[b]).wi, pos.y, precision),
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

"""Given a position, find the most to the left available position without 
overlapping a placed stack."""
function totheleft(pos, solution; precision=3)
    # Find all stacks overlapping on y axis and with x < pos.x
    boxesleft = findboxesleft(pos, Dim(10.0^-precision, 10.0^-precision), solution, precision)

    # Find the stack which extends the most to the right
    rightsides = [get_pos(solution[k]).x + get_dim(solution[k]).le for k in boxesleft]
    leftbound = isempty(rightsides) ? 0 : max(rightsides...)    

    # return the Pos with x position as the right side of the stack
    return Pos(leftbound, pos.y)
end


"""Given a position, find the most to the bottom available position without 
overlapping a placed stack."""
function tothebottom(pos, solution; precision=3)
    # Find all stacks overlapping on x axis and with y < pos.y
    boxesbot = findboxesbelow(pos, Dim(10.0^-precision, 10.0^-precision), solution; precision)

    # Find the stack which extends the most to the top
    topsides = [get_pos(solution[k]).y + get_dim(solution[k]).wi for k in boxesbot]
    botbound = isempty(topsides) ? 0 : max(topsides...)    

    # return the Pos with y position as the top side of the stack
    return Pos(pos.x, botbound)
end

"""
    coveredcorners(corners, o, le, wi; precision=3, verbose=false)

Remove corners covered by provided stack at position o.
"""
function coveredcorners(corners::Vector{<:AbstractPos}, o, le, wi; precision=3, verbose=false)
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

function is_secure(stack, solution; precision=3)
    if eqtol(get_pos(stack).x, 0, precision) 
        return true
    end
    boxesleft = findboxesleft(get_pos(stack), get_dim(stack), solution, precision)
    for k in boxesleft
        b = solution[k]
        if eqtol(get_pos(b).x + get_dim(b).le, get_pos(stack).x)
            return true
        end
    end
    return false
end

function can_be_placed(solution, o::Pos, s::Stack, W, orientation::Symbol; precision=3, verbose=false)

    res = nothing

    if orientation == :Perpendicular
        res = leqtol(o.y + get_dim(s).le, W, precision) && !collision(Pos(o.x, o.y), Dim(get_dim(s).wi, get_dim(s).le), solution; precision)
    elseif orientation == :Parallel
        res = leqtol(o.y + get_dim(s).wi, W, precision) && !collision(Pos(o.x, o.y), Dim(get_dim(s).le, get_dim(s).wi), solution; precision)
    else
        throw(ArgumentError("orientation argument should either be :Perpendicular or :Parallel.\nUnkown flag: $orientation"))
    end 
    return res
    
end

error("must update totheleft and totheright so that they return ProjectedPos corners")

"""
    placestack!(solution::Dict{T, S}, W, i, s::AbstractStack, corners; precision=3, verbose=false) where {T <: Integer, S <: AbstractStack}

Place stack `s` in `solution` in the first available corner in `corners`.
Placing a stack leads to the creation of 2 new corners added to `toadd`.
The corner taken is put in a list `torem` of corners to remove.
"""
function placestack!(solution::Dict{T, S}, W, i, s::AbstractStack, corners::Vector{<:AbstractPos}; precision=3, verbose=false, loading_order=false) where {T <: Integer, S <: AbstractStack}
    torem = []
    toadd = []
    placed = false
    if verbose
        println("Placing s=$s")
    end
    # For each corner potentially available
    for o in corners
        if verbose
            println("Considering corner $o")
        end

        # ignore corners waiting to be removed
        if o in torem
            continue
        end

        orientation = :None

        # if the stack fits oriented with its length perpendicular to the 
        # width of the truck and it doesn't overlap with another stack
        if can_be_placed(solution, get_pos(o), s, W, :Perpendicular; precision)
        # if leqtol(o.y + get_dim(s).le, W, precision) && !collision(Pos(o.x, o.y), Dim(get_dim(s).wi, get_dim(s).le), solution; precision)
            # the stack can be placed in this orientation
            orientation = :Perpendicular
            if verbose
                println("$s can be placed perpendicular")
            end
        # else if the stack fits oriented with its length parallel to the
        # length of the truck and it doesn't overlap with another stack
        elseif can_be_placed(solution, get_pos(o), s, W, :Parallel; precision)
            orientation = :Parallel
            if verbose
                println("$s can be placed parallel")
            end

        else
            if verbose
                println("$s can't be placed at $o because...")
                println("Perpendicular:")
                println("\tleqtol(get_pos(o).y + get_dim(s).le, W, precision) = ", leqtol(get_pos(o).y + get_dim(s).le, W, precision))
                println("\t!collision(Pos(get_pos(o).x, get_pos(o).y), Dim(get_dim(s).wi, get_dim(s).le), solution; precision) = ", !collision(Pos(get_pos(o).x, get_pos(o).y), Dim(get_dim(s).wi, get_dim(s).le), solution; precision, verbose))
                println("Parallel:")
                println("\tleqtol(get_pos(o).y + get_dim(s).wi, W, precision)  = ", leqtol(get_pos(o).y + get_dim(s).wi, W, precision) )
                println("\t!collision(Pos(get_pos(o).x, get_pos(o).y), Dim(get_dim(s).le, get_dim(s).wi), solution; precision) = ", !collision(Pos(get_pos(o).x, get_pos(o).y), Dim(get_dim(s).le, get_dim(s).wi), solution; precision, verbose))

                display(solution)

            end
        end

        if orientation != :None
            placed = true
            if orientation == :Perpendicular
                # Add the stack to this corner in solution
                solution[i] = loading_order ? 
                                                OrderedStack(   Stack(Pos(get_pos(o).x, get_pos(o).y), Dim(get_dim(s).wi, get_dim(s).le)), 
                                                                s.supplier_order, 
                                                                s.supplier_dock_order, 
                                                                s.plant_dock_order) : 
                                                Stack(Pos(get_pos(o).x, get_pos(o).y), Dim(get_dim(s).wi, get_dim(s).le))
                
                # Add new corners
                # corners must be placed as much to the left as possible
                # or as much to the bottom as possible
                push!(toadd,    totheleft(Pos(get_pos(o).x, get_pos(o).y + get_dim(s).le), solution), 
                                tothebottom(Pos(get_pos(o).x + get_dim(s).wi, get_pos(o).y), solution))
            end

            if orientation == :Parallel
                # add to solution
                solution[i] = loading_order ? 
                                                OrderedStack(   Stack(Pos(get_pos(o).x, get_pos(o).y), Dim(get_dim(s).le, get_dim(s).wi)),
                                                                s.supplier_order, 
                                                                s.supplier_dock_order, 
                                                                s.plant_dock_order) : 
                                                Stack(Pos(get_pos(o).x, get_pos(o).y), Dim(get_dim(s).le, get_dim(s).wi))

                # add new corners
                push!(toadd, 
                                totheleft(Pos(get_pos(o).x, get_pos(o).y + get_dim(s).wi), solution), 
                                tothebottom(Pos(get_pos(o).x + get_dim(s).le, get_pos(o).y), solution))
            end

            # remove corner
            push!(torem, o)
            
            # remove covered corners
            covered = coveredcorners(corners, get_pos(o), get_dim(solution[i]).le, get_dim(solution[i]).wi; precision)
            # for o2 in corners
            #     if leqtol(get_pos(o).x, o2.x, precision) && lessertol(o2.x, get_pos(o).x + get_dim(s).le, precision) && leqtol(get_pos(o).y, o2.y, precision) && lessertol(o2.x, get_pos(o).y + get_dim(s).wi, precision)
            #         push!(torem, o2)
            #     end
            
            # end
            
            # update projected corners that are intersected by the newly placed stack
            for o in corners
                if greatertol(get_pos(o).x, get_pos(solution[i]).x + get_dim(solution[i]).le)
                    break
                end
                if is_projected(o) && is_intersected(o, solution[i])
                    push!(to_upd, o)
                end
            end
            for o in to_upd
                upd!(corners, o, solution[i])
            end
            push!(torem, covered...)
            ## It is not enough! Consider the following example:
            """
            ____________________________________
            +---------+
            |         |
            |   a     |
            +---------+
                      :
            +---------X----------+
            |   b     :          |
            +---------:----------+
                      :
            __________O________________________

            b is placed after a. It does not cover the O projected corner.
            Hence, the X projected corner is never added.
            """

            # stop iterating over corners
            break
            
        end
    end

    if !placed
        println("solution\n", solution)
        println("W\n", W)
        println("i\n", i)
        println("s\n", s)
        println("corners\n", corners)
        println("precision\n", precision)
        println("verbose\n", verbose)
        error("The stack can't be placed")
    end

    return torem, toadd
end

"""
    BLtruck(instance::Vector{Pair{T, S}}, W; precision=3, verbose=false, loading_order=false) where {T <: Integer, S <: AbstractStack}

Places stacks in a space of width `W` as to minimize to overall length of the solution.
"""
function BLtruck(instance::Vector{Pair{T, S}}, W; precision=3, verbose=false, loading_order=false) where {T <: Integer, S <: AbstractStack}
    """Lengths must be greater than widths"""
    # TODO pretreatment?
    """One of the two dimensions must be lesser than W?"""
    # TODO

    ## If loading orders must be considered
    if loading_order
        # Sort the stacks
        sort!(instance, by=p -> (supplier_order(p[2]), supplier_dock_order(p[2]), plant_dock_order(p[2])))
    end

    corners = [Pos(0, 0)]
    solution = Dict{Integer, AbstractStack}()
    # torem = Pos[]
    # toadd = Pos[]


    # For each stack to place
    for (i, s) in instance

        if !issorted(corners, by= o -> (get_pos(o).x, get_pos(o).y))
            display(corners)
            error("Not sorted")
        end
        torem, toadd = placestack!(solution, W, i, s, corners; precision=precision, loading_order=loading_order)
        if verbose
            println("About to place stack n. $i, $s")
            println("\tAvailable corners: $corners")
            println("\t$i : $(solution[i]) was placed and added the new corners: $toadd")
            println("\tIt covered the following corners: $torem")
        end
        
        # remove corners waiting to be removed
        filter!(x -> !(x in torem), corners)

        # Add corners to the list of available corners
        push!(corners, toadd...)
        sort!(corners, by= o -> (get_pos(o).x, get_pos(o).y) ) # TODO instead of sorting put elements in right place
        # toadd = Pos[]
        # torem = Pos[]
    end

    return solution
end


"""Return wether position o is illegal with already placed boxes."""
function illegalpos(o, r; precision=3)
    return collision(o, Dim(1/(10^precision), 1/(10^precision)), r)
end

"""
shareside(a::Stack, b::Stack)

Return wether two stacks share a stack.

# Example
```
+---+
| d |
+---+---+---+
| c | a | b |
+-+-+-+-+---+
  | e |
  +---+
```
`a` and `b` share a side.
However, `a` and `e` don't.
"""
function shareside(a::Stack, b::Stack)
    """
    ```
        +---+
        | 1 |
    +---+---+---+
    | 3 | a | 2 |
    +---+---+---+
        | 4 |
        +---+
    ```
    """
    if a.dim.le == b.dim.le
        if  Pos(a.pos.x, a.pos.y + a.dim.wi) == b.pos
            return true
        end

        if Pos(b.pos.x, b.pos.y + b.dim.wi) == a.pos
            return true
        end
    end

    if a.dim.wi == b.dim.wi
        if b.pos == Pos(a.pos.x + a.dim.le, a.pos.y)
            return true
        end
        
        if a.pos == Pos(b.pos.x + b.dim.le, b.pos.y)
            return true
        end
    end

    return false
end

"""
cutrectangle(orientation::String, rectangle::Stack, cut::T; precision=3) where T <: Real

Given a rectangle and depending on orientation reduce the length or width 
in function of the cut, as if the rectangle were cut in two.
Return the new rectangle.
"""
function cutrectangle(orientation::String, rectangle::Stack, cut::T; precision=3) where T <: Real
    res = nothing
    if orientation == "x"
        if eqtol(cut, rectangle.pos.x, precision)
            return rectangle 
        end
        if lessertol(cut, rectangle.pos.x, precision) || geqtol(cut, rectangle.pos.x + rectangle.dim.le, precision)
            throw(ArgumentError("cut ($cut) out of bounds of rectangle $rectangle"))
        end
        res = Stack(rectangle.pos, Dim(cut - rectangle.pos.x , rectangle.dim.wi))
    elseif orientation == "y"
        if eqtol(cut, rectangle.pos.y, precision)
            return rectangle 
        end
        if lessertol(cut, rectangle.pos.y, precision) || geqtol(cut, rectangle.pos.y + rectangle.dim.wi, precision)
            throw(ArgumentError("cut ($cut) out of bounds of rectangle $rectangle"))
        end
        res = Stack(rectangle.pos, Dim(rectangle.dim.le, cut - rectangle.pos.y))
    else
        throw(ArgumentError("orientation $orientation not recognized.\nShould be \"x\" or \"y\"."))
    end
    return res
end

"""
    newrectangle(orientation::String, rectangle::Stack, olddim::Dim, cut::T) where T <: Real

Return the new rectangle obtained by cutting the `rectangle` with the given 
orientation at `cut`.
"""
function newrectangle(orientation::String, rectangle::Stack, olddim::Dim, cut::T) where T <: Real
    if orientation == "x"
        return Stack(Pos(cut, rectangle.pos.y), Dim(rectangle.pos.x + olddim.le - cut, olddim.wi))
    else
        return Stack(Pos(rectangle.pos.x, cut), Dim(olddim.le, rectangle.pos.y + olddim.wi - cut))
    end
end

"""
fuse!(rectangles::Dict{T, Stack}, a::Stack, b::Stack, id) where T <: Integer

Create a new rectangle obtained by fusing the two given rectangles.
Remove old rectangles.
"""
function fuse!(rectangles::Dict{T, Stack}, a::Stack, b::Stack, id) where T <: Integer

    filter!(p -> !(p[2] in [a, b]), rectangles)

    samex = b.pos.x == a.pos.x

    rectangles[id] = Stack(
                                Pos(min(a.pos.x, b.pos.x), min(a.pos.y, b.pos.y)), 
                                samex ? Dim(a.dim.le, a.dim.wi + b.dim.wi) : Dim(a.dim.le + b.dim.le, a.dim.wi))
end

"""
    cutandfuse_generator(L, W, cutiter, fuseiter, precision=3)

This algorithm works in two phases:
1. Have a bigger rectangle, and cut it in two new rectangles.
2. Do the step above a number of time for different rectangles
The two first steps alone are unable to generate some configurations, for instance:
```
+---+-----------+
|   |           |
|   +-------+---+
|   |       |   |
|   |       |   |
+---+-------+   |
|           |   |
+-----------+---+
```
3. Fuse two adjacent compatible (which create a new rectangle) rectangles
4. repeat step 3. a number of times

Example to create the configuration above:
```
+---+-----------+
|   |           |
|   |           |
|   |           |
|   |           |
|   |           |
|   |           |
+---+-----------+

+---+-------+---+
|   |       |   |
|   |       |   |
|   |       |   |
|   |       |   |
|   |       |   |
|   |       |   |
+---+-------+---+

+---+-------+---+
|   |       |   |
+---+-------+---+
|   |       |   |
|   |       |   |
|   |       |   |
|   |       |   |
+---+-------+---+

+---+-------+---+
|   |       |   |
+---+-------+---+
|   |       |   |
|   |       |   |
+---+-------+---+
|   |       |   |
+---+-------+---+
```
And then fuse the right rectangles.

"""
function cutandfuse_generator(L, W, cutiter, fuseiter; precision::Integer=3)
    # println("========")
    rectangles = Dict(0 => Stack(Pos(0, 0), Dim(L, W)))
    nb_rects = 1
    movingL = [L/2]
    movingW = [W/2]

    # cut phase
    for i in 1:cutiter
        # println("---")
        # choose random cut orientation
        orientation = rand(["x", "y"])
        cut = 0.0

        # choose a random cut coordinate
        if orientation == "x"
            cut = rand() * L
            # cut = rand(1:L) # DEBUG

            # The following works
            # cut = rand(movingL)
            # push!(movingL, 1.5 * cut, 0.5 * cut)
            # filter!(x -> x != cut, movingL)
        else
            cut = rand() * W
            # cut = rand(1:W) # DEBUG

            # The following works
            # cut = rand(movingW)
            # push!(movingW, 1.5 * cut, 0.5 * cut)
            # filter!(x -> x != cut, movingW)
        end


        # and find impacted rectangles
        impacted = orientation == "x" ? findboxesabove(Pos(cut, 0), rectangles, precision=precision) :
                                        findboxesright(Pos(0, cut), rectangles, precision=precision)
        
        # display(rectangles)
        # println("cut: $cut on $orientation")
        # println("impacted: $impacted")

        # For every impacted rectangle:
        for k in impacted
            # println("-*-")
            # Update rectangle
            """
            ```
            rectangles[k].pos.x
            __^___
                    rectangles[k].dim.le
                  ____^______
                  +----/----+
                  |    /    |
                  +----/----+
                       ^______ cut
            0 1 2 3 4 5 6 7 8 9 10
            ```
            """
            olddim = Dim(rectangles[k].dim.le, rectangles[k].dim.wi)
            rectangles[k] = cutrectangle(orientation, rectangles[k], cut)
            # println("rectangle $k becomes $(rectangles[k])")
            # if orientation == "x"
            #     # rectangles[k].dim.le = cut - rectangles[k].pos.x 
            #     # rectangles[k].dim = Dim(cut - rectangles[k].pos.x , rectangles[k].dim.wi)
            #     rectangles[k] = Stack(rectangles[k].pos, Dim(cut - rectangles[k].pos.x , rectangles[k].dim.wi))
            # else
            #     # rectangles[k].dim.wi = cut - rectangles[k].pos.y 
            #     # rectangles[k].dim = Dim(rectangles[k].dim.le, cut - rectangles[k].pos.y)
            #     rectangles[k] = Stack(rectangles[k].pos, Dim(rectangles[k].dim.le, cut - rectangles[k].pos.y))
            # end

            # create new formed rectangle only if the cut actually cut the rectangle
            if rectangles[k].dim != olddim
                nb_rects += 1
                rectangles[nb_rects] = newrectangle(orientation, rectangles[k], olddim, cut)
            end
            # println("New rectangle: $(rectangles[nb_rects])")

            # if orientation == "x"
            #     rectangles[i] = Stack(Pos(cut, rectangles[k].pos.y), Dim(rectangles[k].pos.x + olddim.le - cut, olddim.wi))
            # else
            #     rectangles[i] = Stack(Pos(rectangles[k].pos.x, cut), Dim(olddim.le, rectangles[k].pos.y + olddim.wi - cut))
            # end
        end
    end

    # fuse phase
    for i in 1:fuseiter
        # choose a rectangle randomly
        rectangle = rectangles[rand(keys(rectangles))]

        # choose a neighbor randomly
        # TODO optimize
        neighbors = []
        """
        ```
            +---+
            | 1 |
        +---+---+---+
        | 3 |   | 2 |
        +---+---+---+
            | 4 |
            +---+
        ```
        """
        for k in keys(rectangles)
            # if  rectangles[k].pos == Pos(rectangle.pos.x, rectangle.pos.y + rectangle.dim.wi) ||
            #     rectangles[k].pos == Pos(rectangle.pos.x + rectangle.dim.le, rectangle.pos.y) ||
            #     rectangles[k].pos == Pos(rectangle.pos.x, rectangle.pos.y + rectangle.dim.wi)
            #     #...
            # end
            if shareside(rectangle, rectangles[k])
                push!(neighbors, rectangles[k])
            end
        end
        # if no neighbor, continue
        if isempty(neighbors)
            continue
        end
        # choose randomly
        neighbor = rand(neighbors)

        # delete the two rectangles, replace with new bigger one 
        nb_rects += 1
        fuse!(rectangles, rectangle, neighbor, nb_rects)

        # filter!(p -> !(p[2] in [neighbor, rectangle]), rectangles)
        # samex = neighbor.pos.x == rectangle.pos.x
        # rectangles[i+cutiter] = Stack(
        #                             Pos(min(rectangle.pos.x, neighbor.pos.x), min(rectangle.pos.y, neighbor.pos.y)), 
        #                             samex ? Dim(rectangle.dim.le, rectangle.dim.wi + neighbor.dim.wi) : Dim(rectangle.dim.le + neighbor.dim.le, rectangle.dim.wi))

    end

    return rectangles
end

"""
    nb_cuts_fuse_avg(nbstacks)

Return an approximation of the number of necessary cuts to obtain `nbstacks` rectangles.
"""
function nb_cuts_fuse_avg(nbstacks)
    # probability of choosing a side: p[side] = 0.5

    # nbx number of x cuts
    # nby number of y cuts

    # for n the number of iteration, the expected outcome is
    # nbx = nby = n/2

    # each iteration, if the chosen side is x, then the number of new stacks is equal to nby and inversely.
    # Which means that each iteration, the number of new stacks is k/2 for k the current iteration

    """
    nbstacks = 1 + sum_{k=1}^{n} p_side * E[nbx_k] + p_side * E[nby_k] 
    nbstacks = 1 + sum_{k=1}^{n} 0.5 * k/2 + 0.5 * k/2
    nbstacks = 1 + sum_{k=1}^{n} k/2 = 1 + n(n+1)/4
    nbstacks - n(n+1)/4 - 1 = 0
    (4 nbstacks - n(n+1) - 4)/4 = 0
    iff 4 nbstacks - n(n+1) - 4 = 0
    iff - n^2 - n + 4 nbstacks - 4 = 0

    delta = 1 - 4 * -1 * (4 nbstacks - 4) = 1 + 16 nbstacks - 16 = 16 nbstacks - 15
    x1 = (1 + sqrt(16 nbstacks - 15))/-2
    x2 = (1 - sqrt(16 nbstacks - 15))/-2
    """

    return max((1 + sqrt(16 * nbstacks - 15))/-2, (1 - sqrt(16 * nbstacks - 15))/-2)

end
