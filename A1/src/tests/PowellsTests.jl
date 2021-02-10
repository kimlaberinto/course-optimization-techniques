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

xs = -2:0.01:2
ys = @. polynomial_up(xs)
plot(xs,ys)

result = PowellsBracketingMethod(polynomial_up, -4, 0.1, 3)
println(result)

result = PowellsBracketingMethod(polynomial_up, 3, 0.1, 3)
println(result)