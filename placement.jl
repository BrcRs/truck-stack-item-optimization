using Random

# 1 => length
# 2 => width

struct Dim
    le
    wi
end

struct Pos
    x
    y
end

function findboxesabove(pos, dim, r)
    return filter(b -> r[b][1].y >= pos.y && r[b][1].x < pos.x + dim.le && r[b][1].x + r[b][2].le > pos.x, keys(r))
end

function findboxesright(pos, dim, r)
    return filter(b -> r[b][1].x >= pos.x && r[b][1].y < pos.y + dim.wi && r[b][1].y + r[b][2].wi > pos.y, keys(r))
end

function collision(pos, dim, r)
    """Return true if considered box overlaps with existing one"""

    # Find boxes with corresponding x
    # tocheck = filter(b -> r[b][1].x < pos.x + dim.le && r[b][1].x + r[b][2].le > pos.x, keys(r))
    tocheck = findboxesabove(pos, dim, r)
    # check each one
    for k in tocheck
        if pos.y + dim.wi > r[k][1].y + 0.001
            return true
        end
    end

    tocheck = findboxesright(pos, dim, r)
    # check each one
    for k in tocheck
        if pos.x + dim.le > r[k][1].x + 0.001
            return true
        end
    end

    return false
end

function outofbound(pos, dim, W)
    return pos.y + dim.wi > W
end

function place(S, W)
    """Lengths must be greater than widths"""
    # TODO pretreatment
    """One of the two dimensions must be lesser than W"""
    # println("Entering place")
    O = [Pos(0, 0)]
    r = Dict{Any, Any}()
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
            if o.y + s.le <= W && !collision(Pos(o.x, o.y), Dim(s.wi, s.le), r)
                push!(torem, o)
                r[i] = Pos(o.x, o.y), Dim(s.wi, s.le)
                push!(toadd, Pos(o.x, o.y + s.le), Pos(o.x + s.wi, o.y))
                # remove covered starting points
                for o2 in O
                    if o.x <= o2.x < o.x + s.wi && o.y <= o2.y < o.y + s.le
                        push!(torem, o2)
                    end
                end
                break
            elseif o.y + s.wi <= W && !collision(Pos(o.x, o.y), Dim(s.le, s.wi), r)
                push!(torem, o)
                r[i] = Pos(o.x, o.y), Dim(s.le, s.wi)
                push!(toadd, Pos(o.x, o.y + s.wi), Pos(o.x + s.le, o.y))
                # remove covered starting points
                for o2 in O
                    if o.x <= o2.x < o.x + s.le && o.y <= o2.y < o.y + s.wi
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
function genS3(W, L, eps)
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
        println()
        print("#")
        for o in O
            print("=")
            # ro = (ywi = o.y + s.wi, yle =  o.y + s.le, xle = o.x + s.le, xwi = o.x + s.wi)
            wi, le = 0.0, 0.0
            above = findboxesabove(o, Dim(0, 0), r)
            right = findboxesright(o, Dim(0, 0), r)
            lowestyabove = isempty(above) ? W : min([r[k][1].y for k in above]...)
            # generate width
            if lowestyabove - o.y < eps
                print(".")
                wi = lowestyabove - o.y
            else
                # Find lowest limit among already placed boxes
                wi = rand() * (lowestyabove - o.y)
                while o.y + wi > lowestyabove || collision(Pos(o.x, o.y), Dim(0.0, wi), r)
                    print("-")
                    print(lowestyabove-o.y)
                    wi = rand() * (lowestyabove - o.y)
                    print(wi)
                    # sleep(1)
                end
            end
            print("|")
            closestxright = isempty(right) ? L : min([r[k][1].x for k in right]...)
            # generate length
            if closestxright - o.x < eps
                print(".")
                le = closestxright - o.x
            else
                le = rand() * (closestxright - o.x)
                println(collision(Pos(o.x, o.y), Dim(0.0, wi), r))
                println(collision(Pos(o.x, o.y), Dim(le, wi), r))
                while o.x + le > closestxright || collision(Pos(o.x, o.y), Dim(le, wi), r)
                    print("+")
                    # print(closestxright - o.x)
                    le = rand() * (closestxright - o.x)
                    # print(":", le)
                    # println(";", o.x + le)
                    # println("In S:", S)
                    # println("W=", W)
                    # println("L=", L)
                    # display(o)
                    # println("($le, $wi)")
                    # sleep(10)
                end
            end
            
            push!(torem, o)
            push!(toadd, Pos(o.x, o.y + wi), Pos(o.x + le, o.y))
            push!(S, Dim(le, wi))
            r[i] = Pos(o.x, o.y), Dim(le, wi)
        end
        filter!(x -> !(x in torem), O)
        push!(O, toadd...)
        toadd = []
    end

    return S

end

# function genS2(W, L, nbsub, nb)
#     S = []

#     # S = [(step, step) for i in 1:((W*L)/(step*step))]
#     S = Matrix{Integer}(undef, nbsub, nbsub)
#     zones = Dict{Any, Any}()
#     m = 0
#     for i in 1:nbsub
#         for j in 1:nbsub
#             m += 1
#             S[i, j] = m
#             zones[m] = hcat([((i-1) * W/nbsub) ; ((j-1) * L/nbsub)],  [(i * W/nbsub); (j * L/nbsub)])
#         end
#     end
#     tobesel = collect(1:m)
#     while !isempty(tobesel)
#         # Select random cell
#         cell = rand(tobesel)
#         # For four sides
#         for dir in [(-1, 0), (1, 0), (0, 1), (0, -1)]
#             # choose random value for length
#             thelength = rand() * nbsub

#             # cover until size limit or previously expanded rectangle is met
#             while 
#         end
#         filter!(x -> x != cell, tobesel)
#     end

    
#     return S, zones
# end

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
        L = max([r[k][1][1] + r[k][2][1] for k in keys(r)]...)
    end
        # println([r[k][1][1] + r[k][2][1] for k in keys(r)])
    A = Matrix{String}(undef, W, L)
    A .= " "

    for k in keys(r)
        A[r[k][1][2]+1:r[k][1][2]+r[k][2][2], r[k][1][1]+1:r[k][1][1]+r[k][2][1]] .= string(k)
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

function eval(maxW, maxL)
    W = rand(1:maxW)
    L = rand(1:maxL)
    S = rndgenS(W, L, rand(1:min(W, L)))
    r = place(S, W)
    println(W, " ", L)
    println(S)
    println(r)
    foundL = max([r[k][1][1] + r[k][2][1] for k in keys(r)]...)
    return foundL/L
end

function evals(nbiter)
    sum = 0
    for i in 1:nbiter
        sum += eval(20, 20)
    end
    return sum/nbiter
end

