using Revise

using Memento
using Printf
using Test

includet("../A1_module/A1Module.jl")
using .A1Module

setlevel!(getlogger("Main"), "debug")
setlevel!(getlogger(A1Module), "debug")

function parabola_up(x)
    return x^2
end

result = SwannsBracketingMethod(parabola_up, -10, 0.1)
println(result)

result = SwannsBracketingMethod(parabola_up, 10, 0.1)
println(result)

result = SwannsBracketingMethod(parabola_up, 0, 0.1)
println(result)