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

function rosenbrock_banana(input; a = 1, b = .1)
    global N_f_eval += 1
    x = input[1]
    y = input[2]
    return (a - x)^2 + 100 * b * (y - x^2)^2
end

function log_banana(x, y)
    return log(rosenbrock_banana([x, y]) + 1)
end

x = 0:0.01:2
y = 0:0.01:2

X = repeat(reshape(x, 1, :), length(y), 1)
Y = repeat(y, 1, length(x))
p1 = contour(x, y, log_banana, fill = false)
plot(p1)
title!("Q1 Line Search Test (log spaced contours)")

d = [1,1]
x_0 = [0.1, 0.1]
desired_interval_size = 0.1
scatter!([x_0[1]], [x_0[2]], label="x0")


N_f_eval = 0
result = Q1LineSearch(rosenbrock_banana, d, x_0, desired_interval_size; linesearch_method = "SwannsBracketingMethod")
scatter!([result[1]], [result[2]], label="Final (Swanns)")
println("Swanns Result is $result, in N_f_eval=$N_f_eval")

N_f_eval = 0
result = Q1LineSearch(rosenbrock_banana, d, x_0, desired_interval_size; linesearch_method = "PowellsBracketingMethod")
scatter!([result[1]], [result[2]], label="Final (Powells)")
println("Powells Result is $result, in N_f_eval=$N_f_eval")

plot!()

