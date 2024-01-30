include("placement.jl")
include("projected_pos.jl")
include("item.jl")






"""
    placestack!(solution::Dict{T, S}, truck::Truck, i, s::AbstractStack, corners::Vector{<:AbstractPos}; precision=3, verbose=false, loading_order=false) where {T <: Integer, S <: AbstractStack}

Place stack `s` in `solution` in the first available corner in `corners`.
Placing a stack leads to the creation of 2 new corners added to `toadd`.
Additional corners are added when two orthogonal projected corners intersect.
The corner taken is put in a list `torem` of corners to remove.

Specifying `loading_order=true` makes sure OrderedStacks are placed and not standard Stacks.
"""
function placestack!(solution::Dict{T, S}, truck::Truck, i, s::AbstractStack, corners::Vector{<:AbstractPos}; precision=3, verbose=false, loading_order=false) where {T <: Integer, S <: AbstractStack}
    torem = AbstractPos[]
    toadd = AbstractPos[]
    W = get_dim(truck).wi
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
        # + check weight if ItemizedStack
        if can_be_placed(solution, get_pos(o), s, truck, :Perpendicular; precision)
        # if leqtol(o.y + get_dim(s).le, W, precision) && !collision(Pos(o.x, o.y), Dim(get_dim(s).wi, get_dim(s).le), solution; precision)
            # the stack can be placed in this orientation
            orientation = :Perpendicular
            if verbose
                println("$s can be placed perpendicular")
            end
        # else if the stack fits oriented with its length parallel to the
        # length of the truck and it doesn't overlap with another stack
        elseif can_be_placed(solution, get_pos(o), s, truck, :Parallel; precision)
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
                # solution[i] = loading_order ? 
                #                                 OrderedStack(   Stack(Pos(get_pos(o).x, get_pos(o).y), Dim(get_dim(s).wi, get_dim(s).le)), 
                #                                                 get_supplier_order(s), 
                #                                                 get_supplier_dock_order(s), 
                #                                                 get_plant_dock_order(s)) : 
                #                                 Stack(Pos(get_pos(o).x, get_pos(o).y), Dim(get_dim(s).wi, get_dim(s).le)) # DONE have a method create new stack or position_stack common to all AbstractStack
                solution[i] = newStack(s, Pos(get_pos(o).x, get_pos(o).y), Dim(get_dim(s).wi, get_dim(s).le))
                # Add new corners
                # corners must be placed as much to the left as possible
                # or as much to the bottom as possible
                push!(toadd,    totheleft(Pos(get_pos(o).x, get_pos(o).y + get_dim(s).le), solution; precision=precision), 
                                tothebottom(Pos(get_pos(o).x + get_dim(s).wi, get_pos(o).y), solution; precision=precision))
            end

            if orientation == :Parallel
                # add to solution
                # solution[i] = loading_order ? 
                #                                 OrderedStack(   Stack(Pos(get_pos(o).x, get_pos(o).y), Dim(get_dim(s).le, get_dim(s).wi)),
                #                                                 get_supplier_order(s), 
                #                                                 get_supplier_dock_order(s), 
                #                                                 get_plant_dock_order(s)) : 
                #                                 Stack(Pos(get_pos(o).x, get_pos(o).y), Dim(get_dim(s).le, get_dim(s).wi))
                solution[i] = newStack(s, Pos(get_pos(o).x, get_pos(o).y), Dim(get_dim(s).le, get_dim(s).wi))

                # add new corners
                push!(toadd, 
                                totheleft(Pos(get_pos(o).x, get_pos(o).y + get_dim(s).wi), solution; precision=precision), 
                                tothebottom(Pos(get_pos(o).x + get_dim(s).le, get_pos(o).y), solution; precision=precision))
            end



            # compute projected corner intersection and add new corners
            allprojected = filter(x -> is_projected(x), corners)
            toaddprojected = filter(x -> is_projected(x), toadd)
            for c in toaddprojected
                upd_intersection!(toadd, c, convert(Vector{ProjectedPos}, allprojected))
            end

            # remove corner UNLESS IT IS A PROJECTED CORNER
            if !is_projected(o)
                push!(torem, o)
            end
            # remove covered corners
            covered = coveredcorners(corners, get_pos(o), get_dim(solution[i]).le, get_dim(solution[i]).wi; precision)


            # for o2 in corners
            #     if leqtol(get_pos(o).x, o2.x, precision) && lessertol(o2.x, get_pos(o).x + get_dim(s).le, precision) && leqtol(get_pos(o).y, o2.y, precision) && lessertol(o2.x, get_pos(o).y + get_dim(s).wi, precision)
            #         push!(torem, o2)
            #     end
            
            # end
            to_upd = []
            # update projected corners that intersect with the newly placed stack
            for o in corners
                # if greatertol(get_pos(o).x, get_pos(solution[i]).x + get_dim(solution[i]).le)
                #     break
                # end
                if is_projected(o) && is_intersected(o, solution[i])
                    push!(to_upd, o)
                end
            end
            for o in to_upd
                upd!(corners, o, solution[i])
            end
            push!(torem, filter(x -> !is_projected(x), covered)...)
            ## Consider the following example:
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
            Hence, the X projected corner remains but must be updated.
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

# """

# TODO remove since not compatible with how we will deal with weight constraints (add items one by one instead of creating stacks beforehand).
# """
# function BLtruck(items, W, H, plant_dock_orders, supplier_orders, supplier_dock_orders; precision=3, verbose=false)
#     stacks = make_stacks(items, plant_dock_orders, supplier_orders, supplier_dock_orders, H)
#     instance = [Pair(i, s) for (i, s) in enumerate(stacks)]
#     return BLtruck(instance, W; precision=precision, verbose=verbose, loading_order=true, items=true)
# end

"""
    BLtruck(instance::Vector{Pair{T, S}}, truck::Truck; precision=3, verbose=false, loading_order=false) where {T <: Integer, S <: AbstractStack}

Places stacks in a space of width `W` as to minimize to overall length of the solution.

If `loading_order=true` is specified, input stacks are sorted by loading order so that the loading orders constraint is satisfied.
By placing the stacks in order of their loading orders, and due to how the algorithm works, the resulting solution satisfies loading orders.
"""
function BLtruck(instance::Vector{Pair{T, S}}, truck::Truck; precision=3, verbose=false, loading_order=false) where {T <: Integer, S <: AbstractStack}
    println("press enter BLtruck for stacks")
    readline()
    """Lengths must be greater than widths"""
    # TODO pretreatment?
    """One of the two dimensions must be lesser than W?"""
    # TODO
    ## If loading orders must be considered
    if loading_order
        # Sort the stacks
        sort!(instance, by=p -> (get_supplier_order(p[2]), get_supplier_dock_order(p[2]), get_plant_dock_order(p[2])))
    end

    corners = AbstractPos[Pos(0, 0)]
    solution = Dict{Integer, AbstractStack}()
    # torem = Pos[]
    # toadd = Pos[]


    # For each stack to place
    for (i, s) in instance

        if !issorted(corners, by= o -> (get_pos(o).x, get_pos(o).y))
            display(corners)
            error("Not sorted")
        end
        torem, toadd = placestack!(solution, truck, i, s, corners; precision=precision, loading_order=loading_order)
        if verbose
            println("About to place stack n. $i, $s")
            if length(corners) > 20
                # find the position chosen and show neighbor positions
                proxy = first(filter(w -> get_pos(solution[i]).x - get_pos(corners[w]).x < 1, 1:length(corners)))
                println("\t20 Available corners: $(corners[max(proxy-10, 1):min(length(corners), proxy+10)])")
            else
                println("\tAvailable corners: $corners")

            end
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




function BLtruck(instance::Vector{Item}, truck::Truck; precision=3, verbose=false) where {T <: Integer, S <: AbstractStack}
    """Lengths must be greater than widths"""
    # TODO pretreatment?
    """One of the two dimensions must be lesser than W?"""
    # TODO
    println("press enter c")
    readline()

    # Sort the stacks
    sort!(instance, by=item -> (
                        get_supplier_order(truck, get_supplier(item)), 
                        get_supplier_dock_order(truck, get_supplier(item), get_supplier_dock(item)), 
                        get_plant_dock_order(truck, get_plant_dock(item))))
    # sort!(instance, by=item -> (
    #                     get_supplier_orders(truck)[get_supplier(item)], 
    #                     get_supplier_dock_orders(truck)[get_supplier_dock(item)], 
    #                     get_plant_dock_orders(truck)[get_plant_dock(item)]))

    corners = AbstractPos[Pos(0, 0)]
    solution = Dict{Integer, ItemizedStack}()
    # torem = Pos[]
    # toadd = Pos[]


    # For each item to place
    for (i, it) in enumerate(instance)
        println("press enter d", i)
        readline()
        if !issorted(corners, by= o -> (get_pos(o).x, get_pos(o).y))
            display(corners)
            error("Not sorted")
        end
        torem, toadd = placeitem!(solution, truck, it, corners; precision=precision)
        
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

function find_neighbors(rectangles, rectangle)
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
    return neighbors

end
"""
    cutandfuse_generator(L, W, cutiter, fuseiter; precision::Integer=3, mk_squares=false, alpha=0.9, gamma=3)::Dict{Integer, Stack}

Generate a random instance of a placing problem (composed of simple Stack objects). 
It works by first dividing the area of the truck randomly on the X and Y axis and 
then fusing a `fuseiter` number of stacks together to reach complex configurations of solutions.

If `mk_squares=true` is specified, cuts will be made as to harmonize the dimensions of the resulting stacks.
As cutting iterations go on, the cuts become less and less random and will target 
areas where fewer cuts have been made until then. `alpha` specifies how fast the 
transition goes. The closer to 1, the slower. `gamma` also affects the speed at 
which the transition goes. `gamma=1` will result in a linear transition, while 
`gamma>1` will result in a transition that will be faster at the beginning and 
slower as the number of iterations get higher.

# Step by step explanation

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
function cutandfuse_generator(L, W, cutiter, fuseiter; precision::Integer=3, mk_squares=false, alpha=0.9, gamma=3)::Dict{Integer, Stack}
    if alpha > 1 || alpha < 0
        throw(ArgumentError("alpha argument should be between 0 and 1."))
    end
    # println("========")
    rectangles = Dict(0 => Stack(Pos(0, 0), Dim(L, W)))
    nb_rects = 1
    movingL = [L/2]
    movingW = [W/2]

    listL = nothing
    listW = nothing
    if mk_squares
        listL = [(i, v) for (i, v) in enumerate(0:L/10:L)]
        listW = [(i, v) for (i, v) in enumerate(0:W/10:W)]
        weightsL = [0 for i in listL]
        weightsW = [0 for i in listW]
        shuffle!(listL)
        shuffle!(listW)
    end
    step = alpha
    
    # cut phase
    for i in 1:cutiter
        # println("---")
        # choose random cut orientation

        
        orientation = nothing
        if mk_squares
            orient_proba = L / (L + W)
            display(orient_proba)
            if rand() <= orient_proba
                orientation = "x"
            else
                orientation = "y"
            end
        else
            orientation = rand(["x", "y"])
        end

        cut = nothing
        # choose a random cut coordinate
        if orientation == "x"
            if mk_squares
                # Select zone with lowest count
                listL = sort(listL; by=k -> weightsL[k[1]])
                lowest_k = listL[convert(Int64, max(1, min(round(randexp()), L)))]
                display(weightsL)
                display(listL)
                display(lowest_k)
                println()
                # draw random cut around location
                cut = rand() * L * step^gamma + (1-step^gamma) * lowest_k[2]

                # narrow down step
                step *= alpha

                # update weights
                weightsL[lowest_k[1]] += 1

            else
                cut = rand() * L
            end
            # cut = rand(1:L) # DEBUG

            # The following works
            # cut = rand(movingL)
            # push!(movingL, 1.5 * cut, 0.5 * cut)
            # filter!(x -> x != cut, movingL)
        else
            if mk_squares # TODO factorize with L
                # Select zone with lowest count
                listW = sort(listW; by=k -> weightsW[k[1]])
                lowest_k = listW[convert(Int64, max(1, min(round(randexp()), W)))]
                display(weightsW)
                display(listW)
                display(lowest_k)
                println()
                # draw random cut around location
                cut = rand() * W * step^gamma + (1-step^gamma) * lowest_k[2]

                # narrow down step
                step *= alpha

                # update weights
                weightsW[lowest_k[1]] += 1

            else
                cut = rand() * W
            end
            # cut = rand(1:W) # DEBUG

            # The following works
            # cut = rand(movingW)
            # push!(movingW, 1.5 * cut, 0.5 * cut)
            # filter!(x -> x != cut, movingW)
        end


        # and find impacted rectangles
        impacted = orientation == "x" ? findboxesabove(Pos(cut, 0), rectangles, precision=precision) :
                                        findboxesright(Pos(0, cut), rectangles, precision=precision)
        


        # For every impacted rectangle:
        for k in impacted
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

            # create new formed rectangle only if the cut actually cut the rectangle
            if rectangles[k].dim != olddim
                nb_rects += 1
                rectangles[nb_rects] = newrectangle(orientation, rectangles[k], olddim, cut)
            end

        end
    end


    ## Intermediate step before fusing: make sure no rectangle has a dimension smaller than 10.0^-precision
    # do this before fusing to make sure we'll have neighbors
    # as long as such a rectangle exists
    thin_rectangle_k = findfirst(r -> get_dim(r).wi < 10.0^-precision || get_dim(r).le < 10.0^-precision, rectangles)
    while !isnothing(thin_rectangle_k)
        thin_rectangle = rectangles[thin_rectangle_k]
        # fuse it with a neighbor
        neighbors = find_neighbors(rectangles, thin_rectangle)

        # if isempty(neighbors)
        #     display(rectangles)
        #     error("Stack $thin_rectangle_k has no neighbors :(")
        # end

        # DO NOT choose randomly
        # choose the neighbor on the side you want to enlarge
        neighbor = nothing
        if get_dim(thin_rectangle).wi < 10.0^-precision
            neighbors = filter(x -> eqtol(get_dim(x).le, get_dim(thin_rectangle).le, precision), neighbors)
            if isempty(neighbors)
                display(rectangles)
                error("Stack $thin_rectangle_k has no neighbors that share its length :(")
                # this happens once in 5 millions, I don't know why, so I am leaving it here
            end
            neighbor = rand(neighbors)
        else
            neighbors = filter(x -> eqtol(get_dim(x).wi, get_dim(thin_rectangle).wi, precision), neighbors)
            if isempty(neighbors)
                display(rectangles)
                error("Stack $thin_rectangle_k has no neighbors that share its width :(")
            end
            neighbor = rand(neighbors)
        end

        # delete the two rectangles, replace with new bigger one 
        nb_rects += 1
        fuse!(rectangles, thin_rectangle, neighbor, nb_rects)

        
        thin_rectangle_k = findfirst(r -> get_dim(r).wi < 10.0^-precision || get_dim(r).le < 10.0^-precision, rectangles)
    end

    # fuse phase
    for i in 1:fuseiter
        # choose a rectangle randomly
        rectangle = rectangles[rand(keys(rectangles))]

        # choose a neighbor randomly
        # TODO optimize
        neighbors = find_neighbors(rectangles, rectangle)

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
