using Revise

using Plots
using Printf

includet("A1_module/A1Module.jl")
using .A1Module

const ROSENBROCK_A = 1.0
const ROSENBROCK_B = 0.1

N_f_eval = 0
function rosenbrock_banana(x, y; a = ROSENBROCK_A, b = ROSENBROCK_B)
    global N_f_eval += 1
    return (a - x)^2 + 100 * b * (y - x^2)^2
end


begin 
    x_plot = -2:0.01:2
    y_plot = -2:0.01:3
    contour(x_plot, y_plot, rosenbrock_banana, 
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

end

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