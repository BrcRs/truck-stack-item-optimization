using Test

include("../src/placement_visualizer.jl")

solution = Dict(
            0 => Stack(Pos(0, 0), Dim(2.32568, 2.72428)),
            4 => Stack(Pos(0, 2.72428), Dim(97.2757, 0.140078)),
            5 => Stack(Pos(97.2757, 0), Dim(97.2757, 7.53424)),
            6 => Stack(Pos(0, 2.86436), Dim(7.53424, 2.72428)),
            2 => Stack(Pos(0, 5.58863), Dim(0.140078, 2.72428)),
            3 => Stack(Pos(194.551, 0), Dim(97.2757, 2.32568))
        )
W = 10
L = 100
plot_placement(W, L, solution)

NBCUTS = 20
NBFUSE = 40

rectangles = cutandfuse_generator(L, W, NBCUTS, NBFUSE; precision=3)
plot_placement(W, L, rectangles, true)
instance = [pair for pair in rectangles]
solution = BLtruck(instance, W, precision=3)
plot_placement(W, L, solution, true)

# some_rects=[
#            rectangle_from_coords(1, 1, 5, 5)
#            rectangle_from_coords(10, 10, 15, 15)
#        ]

# plot(some_rects[:,1], some_rects[:,2], label = "some group")

# plot(0:5,0:5)
# plot!(rectangle(3,2,0,0), opacity=.5)