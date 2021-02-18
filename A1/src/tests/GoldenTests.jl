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


result, history = GoldenSectionSearch(parabola_up, -10, 10, 0.1)
println(result)

l = @layout [a{0.6h}; b]

xs = -10:0.01:10
ys = @. parabola_up(xs)
p1 = plot(xs,ys)

p2 = plot()
for (i, interval) in enumerate(history)
    plot!([interval[1], interval[2]], [-i, -i], label="$i", shape=:circle, markersize=4, ytick=[], legend=false)
end

plot(p1, p2, layout=l)