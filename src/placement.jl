using Random

# 1 => length
# 2 => width

struct Dim
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
struct Pos
    x
    y
end

struct Stack
    pos::Pos
    dim::Dim
end

function findboxesabove(pos, r, precision=3)
    # return filter(b -> r[b][:pos].y >= pos.y && r[b][:pos].x < pos.x + dim.le && r[b][:pos].x + r[b][:dim].le > pos.x, keys(r))
    # return filter(b -> geqtol(r[b][:pos].y, pos.y, precision) && lessertol(r[b][:pos].x, pos.x + dim.le, precision) && greatertol(r[b][:pos].x + r[b][:dim].le, pos.x, precision), keys(r))
    return filter(b -> geqtol(r[b][:pos].y, pos.y, precision) && leqtol(r[b][:pos].x, pos.x, precision) && greatertol(r[b][:pos].x + r[b][:dim].le, pos.x, precision), keys(r))
end

function findboxesright(pos, r, precision=3)
    # return filter(b -> geqtol(r[b][:pos].x, pos.x, precision) && lessertol(r[b][:pos].y, pos.y + dim.wi, precision) && greatertol(r[b][:pos].y + r[b][:dim].wi, pos.y, precision), keys(r))
    return filter(
        b -> geqtol(r[b][:pos].x, pos.x, precision) && 
        leqtol(r[b][:pos].y, pos.y, precision) && 
        greatertol(r[b][:pos].y + r[b][:dim].wi, pos.y, precision), 
        keys(r))
end

function findboxesright(pos, dim, r, precision=3)
    return filter(
        b -> geqtol(r[b][:pos].x, pos.x, precision) && 
        lessertol(r[b][:pos].y, pos.y + dim.wi, precision) && 
        greatertol(r[b][:pos].y + r[b][:dim].wi, pos.y, precision), 
        keys(r))
    # return filter(
    #     b -> geqtol(r[b][:pos].x, pos.x, precision) && 
    #     leqtol(r[b][:pos].y, pos.y, precision) && 
    #     greatertol(r[b][:pos].y + r[b][:dim].wi, pos.y, precision), 
    #     keys(r))
end

function findboxesleft(pos, dim, r, precision=3)
    return filter(
        b -> leqtol(r[b][:pos].x, pos.x, precision) && 
        lessertol(r[b][:pos].y, pos.y + dim.wi, precision) && 
        greatertol(r[b][:pos].y + r[b][:dim].wi, pos.y, precision),
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
    # tocheck = filter(b -> r[b][:pos].x < pos.x + dim.le && r[b][:pos].x + r[b][:dim].le > pos.x, keys(r))
    # tocheck = findboxesabove(pos, dim, r, precision)
    # check each one
    for k in keys(r)

        if overlapX(pos, dim, r[k][:pos], r[k][:dim]; precision) && overlapY(pos, dim, r[k][:pos], r[k][:dim]; precision)
            if verbose
                println(k)
                println("overlapX($pos, $dim, $(r[k][:pos]), $(r[k][:dim]); $precision) = ", overlapX(pos, dim, r[k][:pos], r[k][:dim]; precision))
                println("overlapY($pos, $dim, $(r[k][:pos]), $(r[k][:dim]); $precision) = ", overlapY(pos, dim, r[k][:pos], r[k][:dim]; precision))
            end
            return true
        end
    end

    return false
end

function outofbound(pos, dim, W, precision=3)
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
            #     println(!collision(Pos(o.x, o.y), Dim(s.wi, s.le), r))
            #     println(!collision(Pos(o.x, o.y), Dim(s.le, s.wi), r))
            # end
            if leqtol(o.y + s.le, W, precision) && !collision(Pos(o.x, o.y), Dim(s.wi, s.le), r; precision)
                push!(torem, o)
                r[i] = Pos(o.x, o.y), Dim(s.wi, s.le)
                push!(toadd, Pos(o.x, o.y + s.le), Pos(o.x + s.wi, o.y))
                # remove covered starting points
                for o2 in O
                    if leqtol(o.x,  o2.x, precision) && lessertol(o2.x, o.x + s.wi, precision) && leqtol(o.y, o2.y, precision) && lessertol(o2.y, o.y + s.le, precision)
                        push!(torem, o2)
                    end
                end
                break
            elseif leqtol(o.y + s.wi, W, precision) && !collision(Pos(o.x, o.y), Dim(s.le, s.wi), r; precision)
                push!(torem, o)
                r[i] = Stack(Pos(o.x, o.y), Dim(s.le, s.wi))
                push!(toadd, Pos(o.x, o.y + s.wi), Pos(o.x + s.le, o.y))
                # remove covered starting points
                for o2 in O
                    if leqtol(o.x, o2.x, precision) && lessertol(o2.x, o.x + s.le, precision) && leqtol(o.y, o2.y, precision) && lessertol(o2.x, o.y + s.wi, precision)
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
    leftbound = max([solution[k][:pos].x + solution[k][:dim].le for k in boxesleft])    

    # return the x position of the right side of the stack
    return leftbound
end

struct Tag
    tag::Symbol
    _tags::Vector{Symbol}
end
function Tag(t)
    _tags = [:Perpendicular, :Parallel, :None]
    if !(t in _tags)
        throw(ArgumentError("$t is not a recognized tag.\nRecognized tags are $_tags"))
    end
    return Tag(t, _tags)
end


function BLtruck(instance, precision=3)
    """Lengths must be greater than widths"""
    # TODO pretreatment?
    """One of the two dimensions must be lesser than W?"""

    corners = [Pos(0, 0)]
    solution = Dict{Any, Any}()
    torem = []
    toadd = []


    # For each stack to place
    for (i, s) in enumerate(instance[:stacks])
    
        # For each corner potentialy available
        for o in corners
    
            # ignore corners waiting to be removed
            if o in torem
                continue
            end

            orientation = Tag(:None)

            # if the stack fits oriented with its length perpendicular to the 
            # width of the truck and it doesn't overlap with another stack
            if leqtol(o.y + s.le, W, precision) && !collision(Pos(o.x, o.y), Dim(s.wi, s.le), solution; precision)
                # the stack can be placed in this orientation
                orientation = Tag(:Perpendicular)
            # else if the stack fits oriented with its length parallel to the
            # length of the truck and it doesn't overlap with another stack
            elseif leqtol(o.y + s.wi, W, precision) && !collision(Pos(o.x, o.y), Dim(s.le, s.wi), solution; precision)
                orientation = Tag(:Parallel)
            end

            if orientation != Tag(:None)
                
                if orientation == Tag(:Perpendicular)
                    # Add the stack to this corner in solution
                    solution[i] = Stack(Pos(o.x, o.y), Dim(s.wi, s.le))
                    
                    # Add new corners
                    # TODO floating corners: corners must be placed as much to the left as possible
                    push!(toadd,    totheleft(Pos(o.x, o.y + s.le), solution), 
                                    totheleft(Pos(o.x + s.wi, o.y), solution))
                end

                if orientation == Tag(:Parallel)
                    # add to solution
                    solution[i] = Stack(Pos(o.x, o.y), Dim(s.le, s.wi))

                    # add new corners
                    # TODO add floating corners
                    push!(toadd, 
                                    totheleft(Pos(o.x, o.y + s.wi), solution), 
                                    totheleft(Pos(o.x + s.le, o.y), solution))
                end

                # remove corner
                push!(torem, o)
                
                # remove covered corners
                for o2 in corners
                    if leqtol(o.x, o2.x, precision) && lessertol(o2.x, o.x + s.le, precision) && leqtol(o.y, o2.y, precision) && lessertol(o2.x, o.y + s.wi, precision)
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

        toadd = []
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
            # ro = (ywi = o.y + s.wi, yle =  o.y + s.le, xle = o.x + s.le, xwi = o.x + s.wi)
            wi, le = missing, missing
            # First draw a vertical line until box is found
            above = findboxesabove(o, r, precision)
            # Find lowest limit among already placed boxes
            lowestyabove = isempty(above) ? W : min([r[k][:pos].y for k in above]...)
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
            closestxright = isempty(right) ? L : min([r[k][:pos].x for k in right]...)
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
            r[i] = Pos(o.x, o.y), Dim(le, wi)
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
        L = max([r[k][:pos][1] + r[k][:dim][1] for k in keys(r)]...)
    end
        # println([r[k][:pos][1] + r[k][:dim][1] for k in keys(r)])
    A = Matrix{String}(undef, W, L)
    A .= " "

    for k in keys(r)
        A[r[k][:pos][2]+1:r[k][:pos][2]+r[k][:dim][2], r[k][:pos][1]+1:r[k][:pos][1]+r[k][:dim][1]] .= string(k)
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
    foundL = max([r[k][:pos].x + r[k][:dim].le for k in keys(r)]...)
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

