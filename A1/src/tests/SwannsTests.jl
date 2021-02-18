using Revise

using Memento
using Printf
using Test
using Plots

includet("../A1_module/A1Module.jl")
using .A1Module

setlevel!(getlogger("Main"), "debug")
setlevel!(getlogger(A1Module), "debug")

function parabola_up(x)
    return x^2
end

function parabola_down(x)
    return -x^2
end


begin
    result, history = SwannsBracketingMethod(parabola_up, -10, 0.1)
    println(result)

    l = @layout [a{0.6h}; b]
    xs = -16:0.01:16
    ys = @. parabola_up(xs)
    p1 = plot(xs,ys)
    xlims!(p1, -16, 16)

    p2 = plot(title="Points A, B, C in each iteration of Swanns")
    for (i, interval) in enumerate(history)
        plot!([interval[1], interval[2], interval[3]], [-i, -i, -i], label="$i", shape=:circle, markersize=4, ytick=[], legend=false)
    end
    xlims!(p2, -16, 16)

    plot(p1, p2, layout=l)
end

begin
    result, history = SwannsBracketingMethod(parabola_up, 10, 0.1)
    println(result)

    l = @layout [a{0.6h}; b]
    xs = -16:0.01:16
    ys = @. parabola_up(xs)
    p1 = plot(xs,ys)
    xlims!(p1, -16, 16)

    p2 = plot(title="Points A, B, C in each iteration of Swanns")
    for (i, interval) in enumerate(history)
        plot!([interval[1], interval[2], interval[3]], [-i, -i, -i], label="$i", shape=:circle, markersize=4, ytick=[], legend=false)
    end
    xlims!(p2, -16, 16)

    plot(p1, p2, layout=l)
end

begin
    result, history = SwannsBracketingMethod(parabola_up, 0, 0.1)
    println(result)

    l = @layout [a{0.6h}; b]
    xs = -2:0.01:2
    ys = @. parabola_up(xs)
    p1 = plot(xs,ys)
    xlims!(p1, -2, 2)

    p2 = plot(title="Points A, B, C in each iteration of Swanns")
    for (i, interval) in enumerate(history)
        plot!([interval[1], interval[2], interval[3]], [-i, -i, -i], label="$i", shape=:circle, markersize=4, ytick=[], legend=false)
    end
    xlims!(p2, -2, 2)

    plot(p1, p2, layout=l)
end

@test_throws ErrorException SwannsBracketingMethod(parabola_down, 0, 0.1)