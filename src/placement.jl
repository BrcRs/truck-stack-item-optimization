using Random
using AutoHashEquals
# 1 => length
# 2 => width

@auto_hash_equals struct Dim
    le # le corresponds to x axis
    wi
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

function findboxesabove(pos, r; precision=3)
    # return filter(b -> r[b].pos.y >= pos.y && r[b].pos.x < pos.x + dim.le && r[b].pos.x + r[b].dim.le > pos.x, keys(r))
    # return filter(b -> geqtol(r[b].pos.y, pos.y, precision) && lessertol(r[b].pos.x, pos.x + dim.le, precision) && greatertol(r[b].pos.x + r[b].dim.le, pos.x, precision), keys(r))
    return filter(b -> geqtol(r[b].pos.y, pos.y, precision) && leqtol(r[b].pos.x, pos.x, precision) && greatertol(r[b].pos.x + r[b].dim.le, pos.x, precision), keys(r))
end

function findboxesright(pos, r; precision::Integer=3)
    # return filter(b -> geqtol(r[b].pos.x, pos.x, precision) && lessertol(r[b].pos.y, pos.y + dim.wi, precision) && greatertol(r[b].pos.y + r[b].dim.wi, pos.y, precision), keys(r))
    return filter(
        b -> geqtol(r[b].pos.x, pos.x, precision) && 
        leqtol(r[b].pos.y, pos.y, precision) && 
        greatertol(r[b].pos.y + r[b].dim.wi, pos.y, precision), 
        keys(r))
end

function findboxesright(pos, dim, r; precision=3)
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

function findboxesleft(pos, dim, r, precision=3)
    return filter(
        b -> leqtol(r[b].pos.x, pos.x, precision) && 
        lessertol(r[b].pos.y, pos.y + dim.wi, precision) && 
        greatertol(r[b].pos.y + r[b].dim.wi, pos.y, precision),
        keys(r))
end

"""Return wether the two boxes overlap on X"""
function overlapX(apos, adim, bpos, bdim; precision=3)
    return greatertol(apos.x + adim.le, bpos.x, precision) && lessertol(apos.x, bpos.x + bdim.le, precision)
end


"""Return wether the two boxes overlap on X"""
function overlapY(apos, adim, bpos, bdim; precision=3)
    return greatertol(apos.y + adim.wi, bpos.y, precision) && lessertol(apos.y, bpos.y + bdim.wi, precision)
end


leqtol(a, b, decimals=3) = round(a, digits=decimals) <= round(b, digits=decimals)

greatertol(a, b, decimals=3) = !eqtol(a, b, decimals) && geqtol(a, b, decimals)
lessertol(a, b, decimals=3) = !eqtol(a, b, decimals) && leqtol(a, b, decimals)

geqtol(a, b, decimals=3) = round(a, digits=decimals) >= round(b, digits=decimals)

eqtol(a, b, decimals=3) = round(a, digits=decimals) == round(b, digits=decimals)

function collision(pos, dim, r; precision=3, verbose=false)
    """Return true if considered box overlaps with existing one"""

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

function place(S, W, precision=3)
    """Lengths must be greater than widths"""
    # TODO pretreatment
    """One of the two dimensions must be lesser than W"""
    # println("Entering place")
    O = [Pos(0, 0)]
    r = Dict{Int, Stack}()
    torem = []
    toadd = []
    for (i, s) in enumerate(S)
        # println(i, " ", s)
        # println("O = ", O)
        for o in O
            if o in torem
                continue
            end
            # if i == 4
            #     println(3 in keys(r))
            #     println(!collision(Pos(o.x, o.y), Dim(s.dim.wi, s.dim.le), r))
            #     println(!collision(Pos(o.x, o.y), Dim(s.dim.le, s.dim.wi), r))
            # end
            if leqtol(o.y + s.dim.le, W, precision) && !collision(Pos(o.x, o.y), Dim(s.dim.wi, s.dim.le), r; precision)
                push!(torem, o)
                r[i] = Stack(Pos(o.x, o.y), Dim(s.dim.wi, s.dim.le))
                push!(toadd, Pos(o.x, o.y + s.dim.le), Pos(o.x + s.dim.wi, o.y))
                # remove covered starting points
                for o2 in O
                    if leqtol(o.x,  o2.x, precision) && lessertol(o2.x, o.x + s.dim.wi, precision) && leqtol(o.y, o2.y, precision) && lessertol(o2.y, o.y + s.dim.le, precision)
                        push!(torem, o2)
                    end
                end
                break
            elseif leqtol(o.y + s.dim.wi, W, precision) && !collision(Pos(o.x, o.y), Dim(s.dim.le, s.dim.wi), r; precision)
                push!(torem, o)
                r[i] = Stack(Pos(o.x, o.y), Dim(s.dim.le, s.dim.wi))
                push!(toadd, Pos(o.x, o.y + s.dim.wi), Pos(o.x + s.dim.le, o.y))
                # remove covered starting points
                for o2 in O
                    if leqtol(o.x, o2.x, precision) && lessertol(o2.x, o.x + s.dim.le, precision) && leqtol(o.y, o2.y, precision) && lessertol(o2.x, o.y + s.dim.wi, precision)
                        push!(torem, o2)
                    end
                end
                break
            end
                
        end
        filter!(x -> !(x in torem), O)
        push!(O, toadd...)
        toadd = []
    end

    return r

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

# struct Tag
#     tag::Symbol
#     _tags::Vector{Symbol}
# end
# function Tag(t)
#     _tags = [:Perpendicular, :Parallel, :None]
#     if !(t in _tags)
#         throw(ArgumentError("$t is not a recognized tag.\nRecognized tags are $_tags"))
#     end
#     return Tag(t, _tags)
# end


function BLtruck(instance::Vector{Pair{T, Stack}}, W; precision=3) where T <: Integer
    """Lengths must be greater than widths"""
    # TODO pretreatment?
    """One of the two dimensions must be lesser than W?"""

    corners = [Pos(0, 0)]
    solution = Dict{Integer, Stack}()
    torem = Pos[]
    toadd = Pos[]


    # For each stack to place
    for (i, s) in instance
    
        # For each corner potentially available
        for o in corners
    
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
            # else if the stack fits oriented with its length parallel to the
            # length of the truck and it doesn't overlap with another stack
            elseif leqtol(o.y + s.dim.wi, W, precision) && !collision(Pos(o.x, o.y), Dim(s.dim.le, s.dim.wi), solution; precision)
                orientation = :Parallel
            end

            if orientation != :None
                
                if orientation == :Perpendicular
                    # Add the stack to this corner in solution
                    solution[i] = Stack(Pos(o.x, o.y), Dim(s.dim.wi, s.dim.le))
                    
                    # Add new corners
                    # TODO floating corners: corners must be placed as much to the left as possible
                    push!(toadd,    totheleft(Pos(o.x, o.y + s.dim.le), solution), 
                                    totheleft(Pos(o.x + s.dim.wi, o.y), solution))
                end

                if orientation == :Parallel
                    # add to solution
                    solution[i] = Stack(Pos(o.x, o.y), Dim(s.dim.le, s.dim.wi))

                    # add new corners
                    # TODO add floating corners
                    push!(toadd, 
                                    totheleft(Pos(o.x, o.y + s.dim.wi), solution), 
                                    totheleft(Pos(o.x + s.dim.le, o.y), solution))
                end

                # remove corner
                push!(torem, o)
                
                # remove covered corners
                for o2 in corners
                    if leqtol(o.x, o2.x, precision) && lessertol(o2.x, o.x + s.dim.le, precision) && leqtol(o.y, o2.y, precision) && lessertol(o2.x, o.y + s.dim.wi, precision)
                        push!(torem, o2)
                    end
                
                end

                # stop iterating over corners
                break
                
            end
        end
        
        # remove corners waiting to be removed
        filter!(x -> !(x in torem), corners)

        # Add corners to the list of available corners
        push!(corners, toadd...)

        toadd = Pos[]
    end

    return solution
end


"""Given a position and a the position of the lowest box above it, return the maximum width possible 
of a box placed at that position."""
function genWidth(o, lowestyabove, eps, precision=3)
    if greatertol(eps, lowestyabove - o.y, precision)
        throw(ArgumentError("Space between $(o.y) and $lowestyabove is less than minimum eps=$eps width."))
    end
    wi = rand() * (lowestyabove - o.y - eps) + eps
    return wi
end

function genLength(o, closestxright, eps; precision=3)
    if greatertol(eps, closestxright - o.x)
        throw(ArgumentError("Space between $(o.x) and $closestxright is less than minimum eps=$eps width."))
    end
    le = rand() * (closestxright - o.x - eps) + eps
    return le
end

"""Return wether position o is illegal with already placed boxes."""
function illegalpos(o, r; precision=3)
    return collision(o, Dim(1/(10^precision), 1/(10^precision)), r)
end

function genS3(W, L, eps, precision=3)
    """Lengths must be greater than widths"""
    """One of the two dimensions must be lesser than W"""
    # println("Entering place")
    O = [Pos(0, 0)]
    S = []
    r = Dict()
    torem = []
    toadd = []
    i = 0
    while !isempty(O)
        i += 1
        # println()
        # print("#")
        for o in O
            if o in torem
                continue
            end
            push!(torem, o)
            if illegalpos(o, r; precision)
                continue
            end
            # print("=")
            # ro = (ywi = o.y + s.dim.wi, yle =  o.y + s.dim.le, xle = o.x + s.dim.le, xwi = o.x + s.dim.wi)
            wi, le = missing, missing
            # First draw a vertical line until box is found
            above = findboxesabove(o, r, precision)
            # Find lowest limit among already placed boxes
            lowestyabove = isempty(above) ? W : min([r[k].pos.y for k in above]...)
            # generate width
            if lessertol(lowestyabove - o.y, eps, precision)
                # print(".")
                wi = lowestyabove - o.y
            else
                wi = genWidth(o, lowestyabove, eps, precision)
            end
            if leqtol(wi, 0, precision)
                continue
            end
            # Sweep to the right to determine closest right box
            right = findboxesright(Pos(o.x, o.y), Dim(precision, wi), r, precision)
            # print("|")
            closestxright = isempty(right) ? L : min([r[k].pos.x for k in right]...)
            # generate length
            if leqtol(closestxright - o.x, eps, precision)
                # print(".")
                le = closestxright - o.x
            else
                le = genLength(o, closestxright, eps; precision)
            end
            
            
            if leqtol(le, 0, precision)
                continue
            end
            # if leqtol(le, 0, precision) || leqtol(wi, 0, precision)
            #     continue
            # end
            push!(toadd, Pos(o.x, o.y + wi), Pos(o.x + le, o.y))
            push!(S, Dim(le, wi))
            r[i] = Stack(Pos(o.x, o.y), Dim(le, wi))
        end
        filter!(x -> !(x in torem), O)
        torem = []
        push!(O, toadd...)
        toadd = []
        # display(r)
        # display(O)
        # sleep(10)
    end

    return S, r

end

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

function cutrectangle(orientation, rectangle, cut)
    res = nothing
    if orientation == "x"
        if cut == rectangle.pos.x
            return rectangle 
        end
        if cut < rectangle.pos.x || cut >= rectangle.pos.x + rectangle.dim.le
            throw(ArgumentError("cut ($cut) out of bounds of rectangle $rectangle"))
        end
        res = Stack(rectangle.pos, Dim(cut - rectangle.pos.x , rectangle.dim.wi))
    else
        if cut == rectangle.pos.y
            return rectangle 
        end
        if cut < rectangle.pos.y || cut >= rectangle.pos.y + rectangle.dim.wi
            throw(ArgumentError("cut ($cut) out of bounds of rectangle $rectangle"))
        end
        res = Stack(rectangle.pos, Dim(rectangle.dim.le, cut - rectangle.pos.y))
    end
    return res
end

function newrectangle(orientation, rectangle, olddim, cut)
    if orientation == "x"
        return Stack(Pos(cut, rectangle.pos.y), Dim(rectangle.pos.x + olddim.le - cut, olddim.wi))
    else
        return Stack(Pos(rectangle.pos.x, cut), Dim(olddim.le, rectangle.pos.y + olddim.wi - cut))
    end
end

function fuse!(rectangles, a, b, id)

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
            # cut = rand() * L
            # cut = rand(1:L) # DEBUG
            cut = rand(movingL)
            push!(movingL, 1.5 * cut, 0.5 * cut)
            filter!(x -> x != cut, movingL)
        else
            # cut = rand() * W
            # cut = rand(1:W) # DEBUG

            cut = rand(movingW)
            push!(movingW, 1.5 * cut, 0.5 * cut)
            filter!(x -> x != cut, movingW)
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

function rndgenS(W, L, ratio)
    println("Entering rndgenS")
    S = []
    pos = []
    r = 0.0
    while r < ratio
        # println(".")
        coord = rand(1:W), rand(1:L)
        dim = rand(1:W), rand(1:L)
        exit = false
        for (i, s) in enumerate(S)
            # println("\t1")
            if (coord[1] >= pos[i][1] && coord[1] <= pos[i][1] + s[1] && coord[2] >= pos[i][2] && coord[1] <= pos[i][2] + s[2]) || pos[i][1] <= coord[1] + dim[1] && pos[i][2] <= coord[2] + dim[2] ||
                coord[1] + dim[1] > W || coord[2] + dim[2] > L
                exit = true
                break
            end
        end
        if exit
            continue
        end
        # status = false
        # dim = (0, 0)
        # while !status
        #     status = true
        #     println("\t2")
        #     print("\t\t")
        #     dim = rand() * W, rand() * L
        #     for (i, s) in enumerate(S)
        #         print(".")
        #         if pos[i][1] <= coord[1] + dim[1] && pos[i][2] <= coord[2] + dim[2] 
        #             status = false
        #             break
        #         end
        #     end
        # end 
        push!(S, dim)
        push!(pos, coord)
        # upd ratio
        r = sum([s[1] * s[2] for s in S]) / (W * L)
        printplacement(Dict{Any, Any}(k => (pos[k], S[k]) for k in eachindex(pos)), W, L)
    end
    return S
end

# W = 2
# ----------


# ----------

function printplacement(r, W, L=nothing)
    if isnothing(L)
        L = max([r[k].pos[1] + r[k].dim[1] for k in keys(r)]...)
    end
        # println([r[k].pos[1] + r[k].dim[1] for k in keys(r)])
    A = Matrix{String}(undef, W, L)
    A .= " "

    for k in keys(r)
        A[r[k].pos[2]+1:r[k].pos[2]+r[k].dim[2], r[k].pos[1]+1:r[k].pos[1]+r[k].dim[1]] .= string(k)
    end

    # display(A)
    println("+", ["-" for u in 1:L]..., "+")
    for i in size(A, 1):-1:1
        print("|")
        for j in 1:size(A, 2)
            print(A[i, j])
        end

        println("|")

    end
    println("+", ["-" for u in 1:L]..., "+")

end

function genS(W, L, n)
    S = []
    X = shuffle(0:L)[1:n]
    Y = shuffle(0:W)[1:n]
    # println(X)
    # println(Y)
    SX = [(i > n ? L : X[i]) - (i-1 >= 1 ? X[i-1] : 0) for i in 1:n+1]
    SY = [(i > n ? W : Y[i]) - (i-1 >= 1 ? Y[i-1] : 0) for i in 1:n+1]
    for sx in SX
        for sy in SY
            push!(S, (sx, sy))
        end

    end
    return filter(x -> x[1] * x[2] > 0, S)
end

function eval(maxW, maxL, shuffleS=false)
    W = rand(1:maxW)
    L = rand(1:maxL)
    S, genr = genS3(W, L, rand(1:min(W, L)))
    if shuffleS
        shuffle!(S)
    end
    r = place(S, W)
    # println(W, " ", L)
    # println(S)
    # println(r)
    foundL = max([r[k].pos.x + r[k].dim.le for k in keys(r)]...)
    if foundL <= L
        println("Singularity")
        println(S)
        println(genr)
        println(r)
        println(L)
        println(foundL)
    end
    return foundL/L
end

function evals(nbiter, shuffleS=false)
    sum = 0
    best = 9999
    worst = 1
    for i in 1:nbiter
        e = eval(20, 20, shuffleS)
        sum += eval(20, 20, shuffleS)
        best = min(best, e)
        worst = max(worst, e)
    end
    return sum/nbiter, best, worst
end

