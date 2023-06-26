include("placement.jl")


function main()
    # W = 2
    # r = place([Dim(1, 1), Dim(1, 1), Dim(2, 1), Dim(1, 2)], W)
    # display(r)
    # printplacement(r, W)
    # W = 3

    # r = place([Dim(1, 1), Dim(1, 1), Dim(3, 3), Dim(2, 1), Dim(1, 1)], W)
    # display(r)
    # begin
    #     for k in keys(r)
    #         println("$k ", !collision(r[k][1], r[k][2], filter(p -> p[1] != k, r)))
    #     end
    # end

    println(evals(100000))
    println(evals(100000, true))
    # S = genS(3, 17, 2)
    # println(S)
    # println(3 * 17)
    # println(sum([s[1] * s[2] for s in S]))
    
    
    # W = 10
    # L = 20
    # res = rndgenS(W, L, 0.3)
    # display(res)
    # r = place(res, W)
    # display(r)
    # printplacement(r, W)


end

main()