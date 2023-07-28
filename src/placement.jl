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
@auto_hash_equals struct Pos
    x
    y
end


struct Stack
    pos::Pos
    dim::Dim
end


"""
    findboxesabove(pos::Pos, r::Dict{T, Stack}; precision=3) where T <: Integer

Return elements of `r` that are directly above and overlapping on x axis with `pos`.
By convention, the right and top borders of a stack aren't counted in the stack.
"""
function findboxesabove(pos::Pos, r::Dict{T, Stack}; precision=3) where T <: Integer
    # return filter(b -> r[b].pos.y >= pos.y && r[b].pos.x < pos.x + dim.le && r[b].pos.x + r[b].dim.le > pos.x, keys(r))
    # return filter(b -> geqtol(r[b].pos.y, pos.y, precision) && lessertol(r[b].pos.x, pos.x + dim.le, precision) && greatertol(r[b].pos.x + r[b].dim.le, pos.x, precision), keys(r))
    return filter(b -> geqtol(r[b].pos.y, pos.y, precision) && leqtol(r[b].pos.x, pos.x, precision) && greatertol(r[b].pos.x + r[b].dim.le, pos.x, precision), keys(r))
end

"""
    findboxesbelow(pos::Pos, dim::Dim, r::Dict{T, Stack}; precision=3) where T <: Integer

Return elements of `r` that are directly below and overlapping on x axis with a 
stack of position `pos` and dimensions `dim`.
By convention, the right and top borders of a stack aren't counted in the stack.
"""
function findboxesbelow(pos::Pos, dim::Dim, r::Dict{T, Stack}; precision=3) where T <: Integer
    # return filter(b -> r[b].pos.y >= pos.y && r[b].pos.x < pos.x + dim.le && r[b].pos.x + r[b].dim.le > pos.x, keys(r))
    # return filter(b -> geqtol(r[b].pos.y, pos.y, precision) && lessertol(r[b].pos.x, pos.x + dim.le, precision) && greatertol(r[b].pos.x + r[b].dim.le, pos.x, precision), keys(r))
    return filter(b ->  leqtol(r[b].pos.y, pos.y, precision) && 
                        leqtol(r[b].pos.x, pos.x + dim.le, precision) && 
                        greatertol(r[b].pos.x + r[b].dim.le, pos.x, precision), keys(r))
end

"""
    findboxesright(pos::Pos, r::Dict{T, Stack}; precision::Integer=3) where T <: Integer

Return elements of `r` that are directly to the right and overlapping on y axis with `pos`.
By convention, the right and top borders of a stack aren't counted in the stack.
"""
function findboxesright(pos::Pos, r::Dict{T, Stack}; precision::Integer=3) where T <: Integer
    # return filter(b -> geqtol(r[b].pos.x, pos.x, precision) && lessertol(r[b].pos.y, pos.y + dim.wi, precision) && greatertol(r[b].pos.y + r[b].dim.wi, pos.y, precision), keys(r))
    return filter(
        b -> geqtol(r[b].pos.x, pos.x, precision) && 
        leqtol(r[b].pos.y, pos.y, precision) && 
        greatertol(r[b].pos.y + r[b].dim.wi, pos.y, precision), 
        keys(r))
end

"""
    findboxesright(pos::Pos, dim::Dim, r::Dict{T, Stack}; precision=3) where T <: Integer

Return elements of `r` that are directly to the right and overlapping on y axis 
with a stack of position `pos` and dimensions `dim`.
By convention, the right and top borders of a stack aren't counted in the stack.
"""
function findboxesright(pos::Pos, dim::Dim, r::Dict{T, Stack}; precision=3) where T <: Integer
    return filter(
        b -> geqtol(r[b].pos.x, pos.x, precision) && 
        lessertol(r[b].pos.y, pos.y + dim.wi, precision) && 
        greatertol(r[b].pos.y + r[b].dim.wi, pos.y, precision), 
        keys(r))
    # return filter(
    #     b -> geqtol(r[b].pos.x, pos.x, precision) && 
    #     leqtol(r[b].pos.y, pos.y, precision) && 
    #     greatertol(r[b].pos.y + r[b].dim.wi, pos.y, precision), 
    #     keys(r))
end

"""
    findboxesleft(pos::Pos, dim::Dim, r::Dict{T, Stack}, precision=3) where T <: Integer

Return elements of `r` that are directly to the left and overlapping on y axis 
with a stack of position `pos` and dimensions `dim`.
By convention, the right and top borders of a stack aren't counted in the stack.
"""
function findboxesleft(pos::Pos, dim::Dim, r::Dict{T, Stack}, precision=3) where T <: Integer
    return filter(
        b -> leqtol(r[b].pos.x, pos.x, precision) && 
        lessertol(r[b].pos.y, pos.y + dim.wi, precision) && 
        greatertol(r[b].pos.y + r[b].dim.wi, pos.y, precision),
        keys(r))
end

"""Return wether two boxes overlap on X"""
function overlapX(apos, adim, bpos, bdim; precision=3)
    return greatertol(apos.x + adim.le, bpos.x, precision) && lessertol(apos.x, bpos.x + bdim.le, precision)
end


"""Return wether two boxes overlap on X"""
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
    collision(pos::Pos, dim::Dim, r::Dict{T, Stack}; precision=3, verbose=false) where T <: Integer

Return true if considered box overlaps with existing one in solution `r`.
"""
function collision(pos::Pos, dim::Dim, r::Dict{T, Stack}; precision=3, verbose=false) where T <: Integer

    # Find boxes with corresponding x
    # tocheck = filter(b -> r[b].pos.x < pos.x + dim.le && r[b].pos.x + r[b].dim.le > pos.x, keys(r))
    # tocheck = findboxesabove(pos, dim, r, precision)
    # check each one
    for k in keys(r)

        if overlapX(pos, dim, r[k].pos, r[k].dim; precision) && overlapY(pos, dim, r[k].pos, r[k].dim; precision)
            if verbose
                println(k)
                println("overlapX($pos, $dim, $(r[k].pos), $(r[k].dim); $precision) = ", overlapX(pos, dim, r[k].pos, r[k].dim; precision))
                println("overlapY($pos, $dim, $(r[k].pos), $(r[k].dim); $precision) = ", overlapY(pos, dim, r[k].pos, r[k].dim; precision))
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
    rightsides = [solution[k].pos.x + solution[k].dim.le for k in boxesleft]
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
    topsides = [solution[k].pos.y + solution[k].dim.wi for k in boxesbot]
    botbound = isempty(topsides) ? 0 : max(topsides...)    

    # return the Pos with y position as the top side of the stack
    return Pos(pos.x, botbound)
end

"""
    coveredcorners(corners, stack; precision=3)

Remove corners covered by provided stack at position o.
"""
function coveredcorners(corners, o, stack; precision=3)
    torem = Pos[]
    for o2 in corners
        if leqtol(o.x, o2.x, precision) && lessertol(o2.x, o.x + stack.dim.le, precision) && leqtol(o.y, o2.y, precision) && lessertol(o2.y, o.y + stack.dim.wi, precision)
            push!(torem, o2)
        end
    
    end
    return torem
end

"""
    placestack!(solution::Dict{T, Stack}, W, i, s::Stack, corners; precision=3, verbose=false) where T <: Integer

Place stack `s` in `solution` in the first available corner in `corners`.
Placing a stack leads to the creation of 2 new corners added to `toadd`.
The corner taken is put in a list `torem` of corners to remove.
"""
function placestack!(solution::Dict{T, Stack}, W, i, s::Stack, corners; precision=3, verbose=false) where T <: Integer
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
        if leqtol(o.y + s.dim.le, W, precision) && !collision(Pos(o.x, o.y), Dim(s.dim.wi, s.dim.le), solution; precision)
            # the stack can be placed in this orientation
            orientation = :Perpendicular
            if verbose
                println("$s can be placed perpendicular")
            end
        # else if the stack fits oriented with its length parallel to the
        # length of the truck and it doesn't overlap with another stack
        elseif leqtol(o.y + s.dim.wi, W, precision) && !collision(Pos(o.x, o.y), Dim(s.dim.le, s.dim.wi), solution; precision)
            orientation = :Parallel
            if verbose
                println("$s can be placed parallel")
            end

        else
            if verbose
                println("$s can't be placed at $o because...")
                println("Perpendicular:")
                println("\tleqtol(o.y + s.dim.le, W, precision) = ", leqtol(o.y + s.dim.le, W, precision))
                println("\t!collision(Pos(o.x, o.y), Dim(s.dim.wi, s.dim.le), solution; precision) = ", !collision(Pos(o.x, o.y), Dim(s.dim.wi, s.dim.le), solution; precision, verbose))
                println("Parallel:")
                println("\tleqtol(o.y + s.dim.wi, W, precision)  = ", leqtol(o.y + s.dim.wi, W, precision) )
                println("\t!collision(Pos(o.x, o.y), Dim(s.dim.le, s.dim.wi), solution; precision) = ", !collision(Pos(o.x, o.y), Dim(s.dim.le, s.dim.wi), solution; precision, verbose))

                display(solution)

            end
        end

        if orientation != :None
            placed = true
            if orientation == :Perpendicular
                # Add the stack to this corner in solution
                solution[i] = Stack(Pos(o.x, o.y), Dim(s.dim.wi, s.dim.le))
                
                # Add new corners
                # corners must be placed as much to the left as possible
                # or as much to the bottom as possible
                push!(toadd,    totheleft(Pos(o.x, o.y + s.dim.le), solution), 
                                tothebottom(Pos(o.x + s.dim.wi, o.y), solution))
            end

            if orientation == :Parallel
                # add to solution
                solution[i] = Stack(Pos(o.x, o.y), Dim(s.dim.le, s.dim.wi))

                # add new corners
                push!(toadd, 
                                totheleft(Pos(o.x, o.y + s.dim.wi), solution), 
                                tothebottom(Pos(o.x + s.dim.le, o.y), solution))
            end

            # remove corner
            push!(torem, o)
            
            # remove covered corners
            push!(torem, coveredcorners(corners, o, s; precision)...)
            # for o2 in corners
            #     if leqtol(o.x, o2.x, precision) && lessertol(o2.x, o.x + s.dim.le, precision) && leqtol(o.y, o2.y, precision) && lessertol(o2.x, o.y + s.dim.wi, precision)
            #         push!(torem, o2)
            #     end
            
            # end

            # stop iterating over corners
            break
            
        end
    end

    if !placed
        error("The stack can't be placed")
    end

    return torem, toadd
end

"""
    BLtruck(instance::Vector{Pair{T, Stack}}, W; precision=3, verbose=false) where T <: Integer

Places stacks in a space of width `W` as to minimize to overall length of the solution.
"""
function BLtruck(instance::Vector{Pair{T, Stack}}, W; precision=3, verbose=false) where T <: Integer
    """Lengths must be greater than widths"""
    # TODO pretreatment?
    """One of the two dimensions must be lesser than W?"""

    corners = [Pos(0, 0)]
    solution = Dict{Integer, Stack}()
    torem = Pos[]
    toadd = Pos[]


    # For each stack to place
    for (i, s) in instance
    
        if !issorted(corners, by= o -> o.x)
            display(corners)
            error("Not sorted")
        end
        torem, toadd = placestack!(solution, W, i, s, corners; precision)
        if verbose
            println("About to place stack n. $i, $s")
            println("Available corners: $corners")
            println("$i : $(solution[i]) was placed and added the new corners: $toadd")
        end
        
        # remove corners waiting to be removed
        filter!(x -> !(x in torem), corners)

        # Add corners to the list of available corners
        push!(corners, toadd...)
        sort!(corners, by= o -> o.x ) # TODO instead of sorting put elements in right place
        toadd = Pos[]
        torem = Pos[]
    end

    return solution
end


# """Given a position and a the position of the lowest box above it, return the maximum width possible 
# of a box placed at that position."""
# function genWidth(o, lowestyabove, eps, precision=3)
#     if greatertol(eps, lowestyabove - o.y, precision)
#         throw(ArgumentError("Space between $(o.y) and $lowestyabove is less than minimum eps=$eps width."))
#     end
#     wi = rand() * (lowestyabove - o.y - eps) + eps
#     return wi
# end

# function genLength(o, closestxright, eps; precision=3)
#     if greatertol(eps, closestxright - o.x)
#         throw(ArgumentError("Space between $(o.x) and $closestxright is less than minimum eps=$eps width."))
#     end
#     le = rand() * (closestxright - o.x - eps) + eps
#     return le
# end

"""Return wether position o is illegal with already placed boxes."""
function illegalpos(o, r; precision=3)
    return collision(o, Dim(1/(10^precision), 1/(10^precision)), r)
end

# function genS3(W, L, eps, precision=3)
#     """Lengths must be greater than widths"""
#     """One of the two dimensions must be lesser than W"""
#     # println("Entering place")
#     O = [Pos(0, 0)]
#     S = []
#     r = Dict()
#     torem = []
#     toadd = []
#     i = 0
#     while !isempty(O)
#         i += 1
#         # println()
#         # print("#")
#         for o in O
#             if o in torem
#                 continue
#             end
#             push!(torem, o)
#             if illegalpos(o, r; precision)
#                 continue
#             end
#             # print("=")
#             # ro = (ywi = o.y + s.dim.wi, yle =  o.y + s.dim.le, xle = o.x + s.dim.le, xwi = o.x + s.dim.wi)
#             wi, le = missing, missing
#             # First draw a vertical line until box is found
#             above = findboxesabove(o, r, precision)
#             # Find lowest limit among already placed boxes
#             lowestyabove = isempty(above) ? W : min([r[k].pos.y for k in above]...)
#             # generate width
#             if lessertol(lowestyabove - o.y, eps, precision)
#                 # print(".")
#                 wi = lowestyabove - o.y
#             else
#                 wi = genWidth(o, lowestyabove, eps, precision)
#             end
#             if leqtol(wi, 0, precision)
#                 continue
#             end
#             # Sweep to the right to determine closest right box
#             right = findboxesright(Pos(o.x, o.y), Dim(precision, wi), r, precision)
#             # print("|")
#             closestxright = isempty(right) ? L : min([r[k].pos.x for k in right]...)
#             # generate length
#             if leqtol(closestxright - o.x, eps, precision)
#                 # print(".")
#                 le = closestxright - o.x
#             else
#                 le = genLength(o, closestxright, eps; precision)
#             end
            
            
#             if leqtol(le, 0, precision)
#                 continue
#             end
#             # if leqtol(le, 0, precision) || leqtol(wi, 0, precision)
#             #     continue
#             # end
#             push!(toadd, Pos(o.x, o.y + wi), Pos(o.x + le, o.y))
#             push!(S, Dim(le, wi))
#             r[i] = Stack(Pos(o.x, o.y), Dim(le, wi))
#         end
#         filter!(x -> !(x in torem), O)
#         torem = []
#         push!(O, toadd...)
#         toadd = []
#         # display(r)
#         # display(O)
#         # sleep(10)
#     end

#     return S, r

# end

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
