using Plots
plotly()
include("placement.jl")

rectangle(x, y, w, h) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
rectangle(stack::AbstractStack) = rectangle(get_pos(stack).x, get_pos(stack).y, get_dim(stack).le, get_dim(stack).wi)

"""
    plot_placement(W, L, solution::Dict{T, S}; orthonormal=false) where {T <: Integer, S <: AbstractStack}

Plot the graph representing a solution of the placement of stacks into a truck 
of width W and length L. The plot will appear in your browser and is interactive 
allowing to zoom, pan, etc. Stack names are displayed in the center of their 
corresponding stacks.
"""
function plot_placement(W, L, solution::Dict{T, S}; orthonormal=false) where {T <: Integer, S <: AbstractStack}
    solution_vector = [s for s in solution]
    rectangles = [rectangle(s) for (i, s) in solution_vector]

    # find furthest x coordinate
    max_x = max([get_pos(s).x + get_dim(s).le for (i, s) in solution_vector]...) * 1.1

    annotations = (

    x=[get_pos(s).x + get_dim(s).le/2 for (i, s) in solution_vector],

    y=[get_pos(s).y + get_dim(s).wi/2 for (i, s) in solution_vector],

    text=[text("Stack $i", 8) for (i, s) in solution_vector],

    )

    # Draw truck boundaries and set main options
    p = plot(rectangle(0, 0, L, W), fillcolor = plot_color(:red, 0.3), 
                    label="Truck", linestyle=:dot, linewidth=5, minorgrid=true, legend=:outertopright, legendfontsize = 7)



    # add rectangles to the plot
    for (i, r) in enumerate(rectangles)
        # plot!(r, color=RGB(rand(0:255)/255, rand(0:255)/255, rand(0:255)/255), label="Stack $i", txt=["Stack $i", "", "", ""])
        plot!(p, r, color=RGB(rand(100:255)/255, rand(100:255)/255, rand(100:255)/255), label="Stack $(solution_vector[i][1])", fontsize=5)
    end
    if orthonormal
        scatter!(p, annotations.x, annotations.y, series_annotations = annotations.text, label="Annotations", ma=0, xlims=(-Inf, max_x), ylims=(-Inf, max_x), tickdir=:out)
    else
        scatter!(p, annotations.x, annotations.y, series_annotations = annotations.text, label="Annotations", ma=0, tickdir=:out)
    end
    # plot!(aspect_ratio=1)
    # display(plot())
    # gui(p)
    display(p)
    # sleep(60)
    # readline()
end
