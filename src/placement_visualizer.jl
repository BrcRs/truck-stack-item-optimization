using Plots

include("placement.jl")

function rectangle_from_coords(xb, yb, xt, yt)
    [
        xb yb
        xt xb
        xt yt
        xb yt
        xb yb
        NaN NaN
    ]
end

function rect_from_stack(s::Stack)
    return rectangle_from_coords(s.pos.x, s.pos.y, s.pos.x + s.dim.le, s.pos.y + s.dim.wi)
end

rectangle(x, y, w, h) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
rectangle(stack::Stack) = rectangle(stack.pos.x, stack.pos.y, stack.dim.le, stack.dim.wi)

function plot_placement(W, L, solution::Dict{T, Stack}) where T <: Integer
    solution_vector = [s for s in solution]
    rectangles = [rectangle(s) for (i, s) in solution_vector]

    trace = (

    x=[s.pos.x + s.dim.le/2 for (i, s) in solution_vector],

    y=[s.pos.y + s.dim.wi/2 for (i, s) in solution_vector],

    text=["Stack $i" for (i, s) in solution_vector],

    )

    plot!(rectangle(0, 0, L, W), fillcolor = plot_color(:red, 0.3), label="Truck", linestyle=:dot, linewidth=5)
    # display(rectangles)
    for (i, r) in enumerate(rectangles)
        # plot!(r, color=RGB(rand(0:255)/255, rand(0:255)/255, rand(0:255)/255), label="Stack $i", txt=["Stack $i", "", "", ""])
        plot!(r, color=RGB(rand(100:255)/255, rand(100:255)/255, rand(100:255)/255), label="Stack $(solution_vector[i][1])", fontsize=5)
    end
    scatter!(trace.x, trace.y, series_annotations = trace.text, label="Annotations", ma=0, )
    # plot!(aspect_ratio=1)
    # display(plot())
    # gui(plot())
end
