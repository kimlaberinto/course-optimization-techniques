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
N_grad_f_eval = 0

function rosenbrock_banana(input; a = ROSENBROCK_A, b = ROSENBROCK_B)
    global N_f_eval += 1
    x = input[1]
    y = input[2]
    return (a - x)^2 + 100 * b * (y - x^2)^2
end

function grad_rosenbrock_banana(input; a = ROSENBROCK_A, b = ROSENBROCK_B)
    global N_grad_f_eval += 1
    x = input[1]
    y = input[2]

    grad_x = -2*a - 400*b*x*(y-x^2) + 2*x # 2*(a - x)*(-1) + (100*b)*(2)*(y-x^2)*(-2*x)
    grad_y = 200*b*(y-x^2) #(100*b)*(2)*(y-x^2)*(1)
    return [grad_x, grad_y]
end

begin 
    x_plot = -2:0.01:2
    y_plot = -2:0.01:3
    contour(x_plot, y_plot, (x, y) -> rosenbrock_banana([x, y]), 
        levels = 0:5:400, 
        fill=false, 
        label = "Rosenbrock",
        c=:black,
        aspect_ratio = :equal,
        size = (600, 600),
        legend = :bottomright)
    
    title!("Example 1D Line Searches on Objective Function")

    plot!([0, 2], [-2, 4], label="A", lw = 3)
    plot!([-2, 4], [2, 0], label="B", lw = 3)
    plot!([-2, 2], [-2, 3], label="C", lw = 3)

    scatter!([1], [1], label="Global Min", shape = :diamond, markersize = 10)

    xlims!(-2, 2)
    ylims!(-2, 3)

    savefig("A1/assets/AllLines.svg")
end

# Function for plot file
function make_bracketing_plot(string_for_plotfile, letter, initial_point, away_point, alpha_minmax, golden_plot_string)
    l = @layout [a{0.45w} [b{0.2h} ; c{0.6h}; d; e]]

    LINE_START = initial_point #Custom
    LINE_DIRECTION = (away_point .- initial_point) #custom
    LINE_DIRECTION = LINE_DIRECTION / norm(LINE_DIRECTION)
    POINT_PLUS_ALPHA = LINE_START .+ 1 .* LINE_DIRECTION
    OneD_function = alpha -> rosenbrock_banana(LINE_START .+ alpha .* LINE_DIRECTION)

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
        scatter!([POINT_PLUS_ALPHA[1]],[POINT_PLUS_ALPHA[2]], label="1 α Away", shape=:circle, markersize = 9)
        println(POINT_PLUS_ALPHA)
        xlims!(-2, 2)
        ylims!(-2, 3)
    end

    #Plot 1D
    begin
        xs = alpha_minmax[1]:0.01:alpha_minmax[2]
        ys = @. OneD_function(xs)
        p_1D = plot(xs,ys, legend=false, title="1D Function")
        xlims!(p_1D, minimum(xs), maximum(xs))
    end

    #Plot 1D Log
    begin
        xs = alpha_minmax[1]:0.01:alpha_minmax[2] #Custom
        ys = @. log10(OneD_function(xs))
        p_1Dlog = plot(xs,ys, legend=false, title = "Log10 1D Function")
        xlims!(p_1Dlog, minimum(xs), maximum(xs))
    end

    # Plot Powell
    begin
        desired_interval_size = 1e-3
        _, powells_bracketing_history, golden_powell_history = Q1LineSearch(rosenbrock_banana, LINE_DIRECTION, LINE_START, desired_interval_size; linesearch_method = "PowellsBracketingMethod")

        num_iterations_powell = length(powells_bracketing_history) 
    
        p_powell = plot(title="$num_iterations_powell Iterations for Powell")
        for (i, interval) in enumerate(powells_bracketing_history)
            plot!([interval[1], interval[2], interval[3]], [-i, -i, -i], label="$i", shape=:circle, markersize=4, ytick=[], legend=false)
        end
        xlims!(p_powell, minimum(xs), maximum(xs))
        ylims!(p_powell, -1*(length(powells_bracketing_history)), 1)
    
    end

    # Plot Swanns
    begin
        _, swanns_bracketing_history, golden_swanns_history  = Q1LineSearch(rosenbrock_banana, LINE_DIRECTION, LINE_START, desired_interval_size; linesearch_method = "SwannsBracketingMethod")
    
        num_iterations_swanns = length(swanns_bracketing_history) 

        p_swanns = plot(title="$num_iterations_swanns Iterations for Swanns")
        for (i, interval) in enumerate(swanns_bracketing_history)
            plot!([interval[1], interval[2], interval[3]], [-i, -i, -i], label="$i", shape=:circle, markersize=4, ytick=[], legend=false)
        end
        xlims!(p_swanns, minimum(xs), maximum(xs))
        ylims!(p_swanns, -1*(length(swanns_bracketing_history)), 1)
    
    end

    xlabel!("Alpha (α)")
    plot(p1, p_1D, p_1Dlog, p_powell, p_swanns, layout=l, size=(800, 500))
    savefig("A1/assets/$string_for_plotfile.svg")

    layout_golden = @layout [a{0.2h} ; [b ; c] [d ; e]]

    #Plot Golden for Powell
    begin
        _, initial_interval = first(golden_powell_history)
        print("Golden Initial = $initial_interval")
        xs = range(initial_interval[1], initial_interval[2], length=100)
        ys = @. log10(OneD_function(xs))
        plot_zoomed_powell  = plot(xs, ys, legend=false, title="Interval found by Powell")
        xlims!(plot_zoomed_powell, minimum(xs), maximum(xs))
    end
    begin
        num_iterations_powell_gold = length(golden_powell_history) 
    
        p_powell_goldensearch = plot(title="$num_iterations_powell_gold Golden Section Search (Powell)")
        for (i, interval) in enumerate(golden_powell_history)
            plot!([interval[1], interval[2], interval[3], interval[4]], [-i, -i, -i, -i], label="$i", shape=:circle, markersize=4, ytick=[], legend=false)
        end
        xlims!(p_powell_goldensearch, minimum(xs), maximum(xs))
        ylims!(p_powell_goldensearch, -1*(length(golden_powell_history)), 1)
    end

    begin
        _, initial_interval = first(golden_swanns_history)
        print("Golden Initial = $initial_interval")
        xs = range(initial_interval[1], initial_interval[2], length=100)
        ys = @. log10(OneD_function(xs))
        plot_zoomed_swanns  = plot(xs, ys, legend=false, title="Interval found by Swanns")
        xlims!(plot_zoomed_swanns, minimum(xs), maximum(xs))
    end
    begin
        num_iterations_powell_gold = length(golden_swanns_history) 
    
        p_swanns_goldensearch = plot(title="$num_iterations_powell_gold Golden Section Search (Swanns)")
        for (i, interval) in enumerate(golden_swanns_history)
            plot!([interval[1], interval[2], interval[3], interval[4]], [-i, -i, -i, -i], label="$i", shape=:circle, markersize=4, ytick=[], legend=false)
        end
        xlims!(p_swanns_goldensearch, minimum(xs), maximum(xs))
        ylims!(p_swanns_goldensearch, -1*(length(golden_swanns_history)), 1)
    end
    title!(p_1Dlog, "Log10 1D Function for Line $letter")
    plot(p_1Dlog, plot_zoomed_powell, p_powell_goldensearch, plot_zoomed_swanns, p_swanns_goldensearch, layout=layout_golden, size=(800, 500))
    savefig("A1/assets/$golden_plot_string.svg")
end

# Block for Q1 Plots
begin
    make_bracketing_plot("LineA_initialbracketing", "A", [0, -2], [2, 4], [-2, 8], "LineA_GoldenComparison")
    make_bracketing_plot("LineB_initialbracketing", "B", [-2, 2], [4, 0], [-2, 8], "LineB_GoldenComparison")
    make_bracketing_plot("LineC_initialbracketing", "C", [-2, -2], [2, 3], [-2, 8], "LineC_GoldenComparison")
end


# Block for Q2 Plots
begin
    x_plot = -2:0.01:2
    y_plot = -2:0.01:3

    plot_Q2_swanns = contour(x_plot, y_plot, (x, y) -> rosenbrock_banana([x, y]), 
        levels = 0:5:400, 
        fill=false, 
        label = "Rosenbrock",
        c=:black,
        legend = :bottomright)
    plot_Q2_powells = contour(x_plot, y_plot, (x, y) -> rosenbrock_banana([x, y]), 
        levels = 0:5:400, 
        fill=false, 
        label = "Rosenbrock",
        c=:black,
        legend = :bottomright)

    plot_Q2_swanns_lowtol = contour(x_plot, y_plot, (x, y) -> rosenbrock_banana([x, y]), 
        levels = 0:5:400, 
        fill=false, 
        label = "Rosenbrock",
        c=:black,
        legend = :bottomright)
    plot_Q2_powells_lowtol = contour(x_plot, y_plot, (x, y) -> rosenbrock_banana([x, y]), 
        levels = 0:5:400, 
        fill=false, 
        label = "Rosenbrock",
        c=:black,
        legend = :bottomright)

    plot_Q2_loss_vs_iter_swanns = plot()
    plot_Q2_loss_vs_iter_powells = plot()
    plot_Q2_loss_vs_iter_swanns_lowtol = plot()
    plot_Q2_loss_vs_iter_powells_lowtol = plot()

    N_f_eval = 0
    N_grad_f_eval = 0
    tolerance = 1e-4
    low_tolerance = 1e-2
    points = [[-2, -2], [-1.5, 1.5], [-1, 3], [-0.5, -1.5], [2, 2]]
    labels = ["D", "E", "F", "G", "H"]
    for (i, (init_point, descent_label)) in enumerate(zip(points, labels))
        result, history = Q2SteepestDescent(rosenbrock_banana, grad_rosenbrock_banana, init_point, tolerance;  linesearch_method = "SwannsBracketingMethod")

        _, all_points = get(history, :Nd_point)
        all_points = transpose(hcat(all_points...))
        plot!(plot_Q2_swanns, all_points[:, 1], all_points[:, 2], label = descent_label, shape=:circle, markersize = 3, lw=3, color=i)
        
        iterations, all_points = get(history, :Nd_point)
        loss = []
        for point2d in all_points
            push!(loss, rosenbrock_banana(point2d))
        end
        plot!(plot_Q2_loss_vs_iter_swanns, iterations, loss, label = descent_label, lw=3, color=i)
    end
    title!(plot_Q2_swanns, "Gradient Descent\nSwanns")
    title!(plot_Q2_loss_vs_iter_swanns, "Obj. Func. Progression\nSwanns")

    for (i, (init_point, descent_label)) in enumerate(zip(points, labels))
        result, history = Q2SteepestDescent(rosenbrock_banana, grad_rosenbrock_banana, init_point, tolerance;  linesearch_method = "PowellsBracketingMethod")

        _, all_points = get(history, :Nd_point)
        all_points = transpose(hcat(all_points...))
        plot!(plot_Q2_powells, all_points[:, 1], all_points[:, 2], label = descent_label, shape=:circle, markersize = 3, lw=3, color=i)
    
        iterations, all_points = get(history, :Nd_point)
        loss = []
        for point2d in all_points
            push!(loss, rosenbrock_banana(point2d))
        end
        plot!(plot_Q2_loss_vs_iter_powells, iterations, loss, label = descent_label, lw=3, color=i)
    end
    title!(plot_Q2_powells, "Gradient Descent\nPowells")
    title!(plot_Q2_loss_vs_iter_powells, "Obj. Func. Progression\nPowell")

    for (i, (init_point, descent_label)) in enumerate(zip(points, labels))
        result, history = Q2SteepestDescent(rosenbrock_banana, grad_rosenbrock_banana, init_point, low_tolerance;  linesearch_method = "SwannsBracketingMethod")

        _, all_points = get(history, :Nd_point)
        all_points = transpose(hcat(all_points...))
        plot!(plot_Q2_swanns_lowtol, all_points[:, 1], all_points[:, 2], label = descent_label, shape=:circle, markersize = 3, lw=3, color=i)
    
        iterations, all_points = get(history, :Nd_point)
        loss = []
        for point2d in all_points
            push!(loss, rosenbrock_banana(point2d))
        end
        plot!(plot_Q2_loss_vs_iter_swanns_lowtol, iterations, loss, label = descent_label, lw=3, color=i)
    end
    title!(plot_Q2_swanns_lowtol, "Gradient Descent (low tol.)\nSwanns")
    title!(plot_Q2_loss_vs_iter_swanns_lowtol, "Obj. Func. Progression\nSwanns (low tol.)")

    for (i, (init_point, descent_label)) in enumerate(zip(points, labels))
        result, history = Q2SteepestDescent(rosenbrock_banana, grad_rosenbrock_banana, init_point, low_tolerance;  linesearch_method = "PowellsBracketingMethod")

        _, all_points = get(history, :Nd_point)
        all_points = transpose(hcat(all_points...))
        plot!(plot_Q2_powells_lowtol, all_points[:, 1], all_points[:, 2], label = descent_label, shape=:circle, markersize = 3, lw=3, color=i)
        
        iterations, all_points = get(history, :Nd_point)
        loss = []
        for point2d in all_points
            push!(loss, rosenbrock_banana(point2d))
        end
        plot!(plot_Q2_loss_vs_iter_powells_lowtol, iterations, loss, label = descent_label, lw=3, color=i)
    end
    title!(plot_Q2_powells_lowtol, "Gradient Descent (low tol.)\nPowells")
    title!(plot_Q2_loss_vs_iter_powells_lowtol, "Obj. Func. Progression\nPowells (low tol.)")

    layout_Q2 =  @layout [a b; c d]
    plot(plot_Q2_swanns, plot_Q2_powells, plot_Q2_swanns_lowtol, plot_Q2_powells_lowtol, layout = layout_Q2, size=(1000, 1000))
    savefig("A1/assets/Q2_stepsvisualized.svg")

    layout_Q2_loss_vs_steps = @layout [a b; c d]
    plot(plot_Q2_loss_vs_iter_swanns, plot_Q2_loss_vs_iter_powells, plot_Q2_loss_vs_iter_swanns_lowtol, plot_Q2_loss_vs_iter_powells_lowtol, layout = layout_Q2_loss_vs_steps, legend=true, size=(650, 650))
    plot!(yscale=:log10)
    xlabel!("Number of Gradient Steps")
    ylabel!("Objective Function Value (log scale)")
    savefig("A1/assets/Q2_loss_vs_steps.svg")
end