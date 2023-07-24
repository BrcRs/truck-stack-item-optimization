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

    # println(collision(Pos(0.992648057, 0.30824860456), Dim(0.001, 0.001), Dict(1 => (Pos(0.992648, 0), Dim(0.00735194, 0.97499)))))


    volume = 0.0
    w = 1
    le = 1
    e = 0.01
    S, r = genS3(w, le, e)
    begin
        for k in keys(r)
            # r = Dict(
            #     i => (s.le, s.wi) for (i, s) in enumerate(S))
            println("k: ", !collision(r[k][1], r[k][2], filter(p -> p[1] != k, r)))
            volume += r[k][2].le * r[k][2].wi
        end
    end
    # testoutofbound(r, w)
    # Filling test
    prinlnt("vol ratio ", volume/(w*le) > 0.9)

    # println(evals(100000))
    # println(evals(100000, true))
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