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

    NBCUTS = nb_cuts_fuse_avg(4)
    ITERS = 100000
    nb = 0
    for i in 1:ITERS
        rectangles = cutandfuse_generator(L, W, NBCUTS, 0; precision=3)
        nb += length(rectangles)
    end
    display(nb/ITERS)
    error()


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



    """
    0 => Stack(Pos(0, 0), Dim(80.2629, 5.0879))
    4 => Stack(Pos(0, 5.0879), Dim(19.7371, 1.30617))
    5 => Stack(Pos(80.2629, 0), Dim(80.2629, 3.60593))
    6 => Stack(Pos(0, 6.39407), Dim(19.7371, 3.60593))
    2 => Stack(Pos(80.2629, 3.60593), Dim(19.7371, 5.0879))
    3 => Stack(Pos(160.526, 0), Dim(80.2629, 1.30617))

    Dict{Integer, Stack} with 6 entries:
    0 => Stack(Pos(0, 0), Dim(80, 5))
    4 => Stack(Pos(0, 5), Dim(20, 1))
    5 => Stack(Pos(80, 0), Dim(80, 4))
    6 => Stack(Pos(0, 6), Dim(20, 4))
    2 => Stack(Pos(80, 4), Dim(20, 5))
    3 => Stack(Pos(161, 0), Dim(80, 1))
    
   |____________________________________________________________________________________________________
   |66666666666666666666
   |6                  6                                                           22222222222222222222
   |6                  6                                                           2                  2
   |66666666666666666666                                                           2                  2
   |44444444444444444444                                                           2                  2
   |000000000000000000000000000000000000000000000000000000000000000000000000000000022222222222222222222 # Why 3 wasn't placed here?
   |0                                                                              55555555555555555555555555555555555555555555555555555555555555555555555555555555
   |0                                                                              5                                                                              5
   |0                                                                              5                                                                              5
   |00000000000000000000000000000000000000000000000000000000000000000000000000000005555555555555555555555555555555555555555555555555555555555555555555555555555555533333333333333333333333333333333333333333333333333333333333333333333333333333333
   |____________________________________________________________________________________________________


   Dict{Integer, Stack} with 6 entries:
  0 => Stack(Pos(0, 0), Dim(2.37662, 2.1253))
  4 => Stack(Pos(0, 2.1253), Dim(97.8747, 0.441835))
  5 => Stack(Pos(0, 2.56714), Dim(7.18154, 2.1253))
  6 => Stack(Pos(7.18154, 2.56714), Dim(97.8747, 7.18154))
  2 => Stack(Pos(0, 4.69244), Dim(0.441835, 2.1253))
  3 => Stack(Pos(105.056, 2.56714), Dim(97.8747, 2.37662))

   Dict{Integer, Stack} with 6 entries:
   0 => Stack(Pos(0, 0), Dim(2, 2))
   4 => Stack(Pos(0, 2), Dim(98, 0))
   5 => Stack(Pos(0, 3), Dim(7, 2))
   6 => Stack(Pos(7, 3), Dim(98, 7))
   2 => Stack(Pos(0, 5), Dim(0, 2))
   3 => Stack(Pos(105, 3), Dim(98, 2))

   |____________________________________________________________________________________________________
   |       66666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666
   |       6                                                                                                6
   |       6                                                                                                6
   |2      6                                                                                                6
   |2      6                                                                                                6
   |55555556                                                                                                633333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
   |55555556666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666633333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
   |44444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
   |00
   |00                                                                                                      
   |____________________________________________________________________________________________________


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