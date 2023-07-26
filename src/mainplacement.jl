using Test
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
    # for i in 1:1000
    #     solution = cutandfuse_generator(10, 10, 4, 4, precision=3)
    #     filling = sum([solution[k].dim.le * solution[k].dim.wi for k in keys(solution)])
    #     if filling != 100
    #         display(solution)
    #         println(filling)
    #         break
    #     end
    # end
    L = 1000
    W = 100
    rectangles = cutandfuse_generator(L, W, 2, 2; precision=3)
    instance = [pair for pair in rectangles]
    solution = BLtruck(instance, W, precision=3)
    display(solution)
    foundL = max([solution[k].pos.x + solution[k].dim.le for k in keys(solution)]...)
    display(foundL)
    @testset "No overlapping" begin
        for (j, stack) in solution
            @test !collision(stack.pos, stack.dim, filter(p -> p[1] != j, solution); precision=3, verbose=false)
        end
    end
    @testset "Not out of bounds" begin
        
        for (j, stack) in solution
            @test !outofbound(stack.pos, stack.dim, W; precision=3)
        end
    end

    """

    +------------------+
    |     6            |
    +------------------+
    |         5        |
    +------------------+

    """


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