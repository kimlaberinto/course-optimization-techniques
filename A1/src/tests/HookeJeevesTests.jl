using Revise

using Memento
using Plots
using Printf
using Test

includet("../A1_module/A1Module.jl")
using .A1Module

setlevel!(getlogger("Main"), "debug")
setlevel!(getlogger(A1Module), "debug")
Memento.config!("debug")

const ROSENBROCK_A = 1
const ROSENBROCK_B = .1

N_f_eval = 0
N_grad_f_eval = 0

function rosenbrock_banana(input; a = ROSENBROCK_A, b = ROSENBROCK_B)
    global N_f_eval += 1
    x = input[1]
    y = input[2]
    return (a - x)^2 + 100 * b * (y - x^2)^2
end

result, history = HookeJeeves(rosenbrock_banana, [-2, -1], .1, 0.000001, [[1, 0], [0, 1]])

x_plot = -2:0.01:2
y_plot = -2:0.01:3
plot_path_RB = contour(x_plot, y_plot, (x, y) -> rosenbrock_banana([x, y]), 
    levels = 0:5:400, 
    fill=false, 
    label = "Rosenbrock",
    c=:black,
    aspect_ratio = :equal,
    size = (600, 600),
    legend = :bottomright)

is, all_points = get(history, :x_1)
all_points = transpose(hcat(all_points...))
plot_path_RB = plot!(plot_path_RB, all_points[:, 1], all_points[:, 2], shape=:circle, markersize = 3, lw=3)

errors = []
for (i, current_point) in enumerate(history, :x_1)
    error = rosenbrock_banana(current_point)
    push!(errors, error)
end

plot_lossvssteps_RB = plot(is, errors, yscale=:log10, lw=3, shape=:circle, markersize = 3)
title!(plot_lossvssteps_RB, "Obj. Func. vs HJ Steps")
xlabel!(plot_lossvssteps_RB, "Hooke Jeeves Outer-loop Iterations")
ylabel!(plot_lossvssteps_RB, "Rosenbrock Banana Function Value")

layout_RB = @layout [a{0.5w} b]
plot(plot_path_RB, plot_lossvssteps_RB, layout=layout_RB, size=(700, 400))