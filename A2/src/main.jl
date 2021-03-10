include("makeplots.jl")

println("\nHookeJeeves")
@time evaluateHookeJeeves()
println("\nGradient Descent")
@time evaluateGradientDescent()