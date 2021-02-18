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

function polynomial_up(x)
    return 0.1*x^4 + 0.5x^3 + 0.6x^2 + (-0.8)x + 0
end


begin
    result, history = PowellsBracketingMethod(polynomial_up, -4, 0.1, 3)
    println(result)

    l = @layout [a{0.6h}; b]
    xs = -5:0.01:2
    ys = @. polynomial_up(xs)
    p1 = plot(xs,ys)
    xlims!(p1, -5, 2)

    p2 = plot(title="Points A, B, C in each iteration of Powell")
    for (i, interval) in enumerate(history)
        plot!([interval[1], interval[2], interval[3]], [-i, -i, -i], label="$i", shape=:circle, markersize=4, ytick=[], legend=false)
    end
    xlims!(p2, -5, 2)

    plot(p1, p2, layout=l)
end


begin
    result, history = PowellsBracketingMethod(polynomial_up, 3, 0.1, 3)
    println(result)


    l = @layout [a{0.6h}; b]
    xs = -1:0.01:4
    ys = @. polynomial_up(xs)
    p1 = plot(xs,ys)
    xlims!(p1, -1, 4)

    p2 = plot(title="Points A, B, C in each iteration of Powell")
    for (i, interval) in enumerate(history)
        plot!([interval[1], interval[2], interval[3]], [-i, -i, -i], label="$i", shape=:circle, markersize=4, ytick=[], legend=false)
    end
    xlims!(p2, -1, 4)

    plot(p1, p2, layout=l)
end