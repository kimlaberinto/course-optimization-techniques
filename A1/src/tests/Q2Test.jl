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

function rosenbrock_banana(x, y; a = ROSENBROCK_A, b = ROSENBROCK_B)
    global N_f_eval += 1
    return (a - x)^2 + 100 * b * (y - x^2)^2
end

function grad_rosenbrock_banana(x, y; a = ROSENBROCK_A, b = ROSENBROCK_B)
    global N_grad_f_eval += 1
    grad_x = -2*a - 400*b*x*(y-x^2) + 2*x # 2*(a - x)*(-1) + (100*b)*(2)*(y-x^2)*(-2*x)
    grad_y = 200*b*(y-x^2) #(100*b)*(2)*(y-x^2)*(1)
    return [grad_x, grad_y]
end

function log_banana(x, y)
    return log(rosenbrock_banana(x, y) + 1)
end

x = -4:0.01:4
y = -1:0.01:4

X = repeat(reshape(x, 1, :), length(y), 1)
Y = repeat(y, 1, length(x))
p1 = contour(x, y, log_banana, fill = false)
plot(p1)
title!("Log spaced contours")

N_f_eval = 0
N_grad_f_eval = 0
result, history = Q2SteepestDescent(rosenbrock_banana, grad_rosenbrock_banana, [0, 0], 1e-4)

for point in history
    scatter!([point[1]], [point[2]])
end

plot!()