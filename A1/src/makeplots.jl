using Revise

using Plots
using Printf
using LinearAlgebra
using LaTeXStrings


includet("A1_module/A1Module.jl")
using .A1Module

const ROSENBROCK_A = 1.0
const ROSENBROCK_B = 0.1

N_f_eval = 0
function rosenbrock_banana(input; a = ROSENBROCK_A, b = ROSENBROCK_B)
    global N_f_eval += 1
    x = input[1]
    y = input[2]
    return (a - x)^2 + 100 * b * (y - x^2)^2
end


begin 
    x_plot = -2:0.01:2
    y_plot = -2:0.01:3
    contour(x_plot, y_plot, (x, y) -> rosenbrock_banana([x, y]), 
        levels = 0:5:400, 
        fill=false, 
        label = "Rosenbrock",
        c=cgrad(:thermal, scale=:exp, rev=false),
        aspect_ratio = :equal,
        size = (600, 600),
        legend = :bottomright)
    
    title!("Example 1D Line Searches on Objective Function")

    plot!([0, 2], [-2, 4], label="A", lw = 3)
    plot!([-2, 4], [2, 0], label="B", lw = 3)
    plot!([-2, 2], [-1, 0], label="C", lw = 3)
    plot!([-2, 0], [0, 3], label="D", lw = 3)

    scatter!([1], [1], label="Global Min", shape = :diamond, markersize = 10)

    xlims!(-2, 2)
    ylims!(-2, 3)

    savefig("A1/assets/AllLines.svg")
end

# Function for plot file
function make_bracketing_plot(string_for_plotfile, letter, initial_point, away_point, alpha_minmax)
    l = @layout [a{0.45w} [b{0.2h} ; c{0.6h}; d; e]]

    LINE_START = initial_point #Custom
    LINE_DIRECTION = (away_point .- initial_point) #custom
    println("$initial_point $LINE_DIRECTION")
    LINE_DIRECTION = LINE_DIRECTION / norm(LINE_DIRECTION)
    println("$initial_point $LINE_DIRECTION")
    B0_PLUS_ALPHA = LINE_START .+ 1 .* LINE_DIRECTION
    OneD_LineB_function = alpha -> rosenbrock_banana(LINE_START .+ alpha .* LINE_DIRECTION)

    begin
        x_plot = -2:0.01:2
        y_plot = -2:0.01:3
        p1 = contour(x_plot, y_plot, (x, y) -> rosenbrock_banana([x, y]), 
            levels = 0:5:400, 
            fill=false, 
            label = "Rosenbrock",
            c=:black,
            legend = :bottomright)

        plot!([initial_point[1], away_point[1]], [initial_point[2], away_point[2]], label=letter, lw = 3, title="Line $letter Initial Bracketing")
        scatter!([LINE_START[1]],[LINE_START[2]], label="$letter Init", shape=:circle, markersize = 9)
        scatter!([B0_PLUS_ALPHA[1]],[B0_PLUS_ALPHA[2]], label="1 α Away", shape=:circle, markersize = 9)
        println(B0_PLUS_ALPHA)
        xlims!(-2, 2)
        ylims!(-2, 3)
    end

    #Plot 1D
    begin
        xs = alpha_minmax[1]:0.01:alpha_minmax[2]
        ys = @. OneD_LineB_function(xs)
        p_1D = plot(xs,ys, legend=false, title="1D Function")
        xlims!(p_1D, minimum(xs), maximum(xs))
    end

    #Plot 1D Log
    begin
        xs = alpha_minmax[1]:0.01:alpha_minmax[2] #Custom
        ys = @. log10(OneD_LineB_function(xs))
        p_1Dlog = plot(xs,ys, legend=false, title = "Log 1D Function")
        xlims!(p_1Dlog, minimum(xs), maximum(xs))
    end

    # Plot Powell
    begin
        result, history = PowellsBracketingMethod(OneD_LineB_function, 0, 1, 16)
        println(result)
    
        p_powell = plot(title="Iterations for Powell")
        for (i, interval) in enumerate(history)
            plot!([interval[1], interval[2], interval[3]], [-i, -i, -i], label="$i", shape=:circle, markersize=4, ytick=[], legend=false)
        end
        xlims!(p_powell, minimum(xs), maximum(xs))
        ylims!(p_powell, -1*(length(history)), 1)
    
    end

    # Plot Swanns
    begin
        result, history = SwannsBracketingMethod(OneD_LineB_function, 0, 1)
        println(result)
    
        p_swanns = plot(title="Iterations for Swanns")
        for (i, interval) in enumerate(history)
            plot!([interval[1], interval[2], interval[3]], [-i, -i, -i], label="$i", shape=:circle, markersize=4, ytick=[], legend=false)
        end
        xlims!(p_swanns, minimum(xs), maximum(xs))
        ylims!(p_swanns, -1*(length(history)), 1)
    
    end

    xlabel!("Alpha (α)")
    plot(p1, p_1D, p_1Dlog, p_powell, p_swanns, layout=l, size=(800, 500))
    savefig("A1/assets/$string_for_plotfile.svg")
end

begin
    make_bracketing_plot("LineB_initialbracketing", "B", [-2, 2], [4, 0], [-2, 8])
end

# # Plot for Line B
# begin
#     l = @layout [ [p1 ; p2] [p3 ; p4] ]
# end

# A way to visualize one run of Q1, need Swann's vs Powell
# Left 
begin
    plot()
    xlabel!("Alpha")
    ylabel!("Objective Function Value")
    title!("1D Bracket Finding - Line A - f(x0 + alpha*d)")
end

# Objective Function vs Iteration (DRAFT)
if false
    is, points = get(history)
    yvals = @. rosenbrock_banana(points)
    plot(is, yvals, yscale=:log10, shape=:circle, markersize=3)
    title!("Objective Function Value vs Gradient Descent Iteration")
    xlabel!("Number of Gradient Descent Iterations")
    ylabel!("Rosenbrock Banana Function Value")
end