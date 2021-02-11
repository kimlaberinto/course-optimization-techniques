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

N_f_eval = 0

function rosenbrock_banana(x, y; a = 1, b = .1)
    global N_f_eval += 1
    return (a - x)^2 + 100 * b * (y - x^2)^2
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

d = [1, 1]
x_0 = [0.75, 0.75]
desired_interval_size = 0.1
N_f_eval = 0
result = Q1LineSearch(rosenbrock_banana, d, x_0, desired_interval_size)

scatter!([x_0[1] result[1]], [x_0[2] result[2]])