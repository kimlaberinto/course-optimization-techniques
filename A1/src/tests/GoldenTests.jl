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

function parabola_up(x)
    return x^2
end

# xs = -2:0.01:2
# ys = @. parabola_up(xs)
# plot(xs,ys)

result = GoldenSectionSearch(parabola_up, -10, 10, 0.1)
println(result)