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

N_f_eval = 1
N_grad_f_eval = 1

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

function log_banana(x, y)
    return log(rosenbrock_banana([x, y]) + 1)
end

x = -3:0.01:3
y = -2:0.01:2

X = repeat(reshape(x, 1, :), length(y), 1)
Y = repeat(y, 1, length(x))
p1 = contour(x, y, log_banana, fill = false)
plot(p1, legend=false, aspect_ratio = :equal)
title!("Q2 Steepest Descent (log spaced contours)")

N_f_eval = 0
N_grad_f_eval = 0
result, history = Q2SteepestDescent(rosenbrock_banana, grad_rosenbrock_banana, [3, -1], 1e-4;  linesearch_method = "SwannsBracketingMethod")


_, all_points = get(history, :Nd_point)
all_points = transpose(hcat(all_points...))
plot!(all_points[:, 1], all_points[:, 2], shape=:circle, markersize = 3, lw=3)

plot!()