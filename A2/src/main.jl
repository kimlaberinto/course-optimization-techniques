include("makeplots.jl")

println("\nHookeJeeves")
@time evaluateHookeJeeves()

println("\nGradient Descent")
@time evaluateGradientDescent()

println("\nConjugate Gradient (Fletcher-Reeves)")
@time evaluateConjugateGradientFletcherReeves()

println("\nConjugate Gradient (Hestenes-Stiefel)")
@time evaluateConjugateGradientHestenesStiefel()

println("\nConjugate Gradient (Polak-Ribiere)")
@time evaluateConjugateGradientPolakRibiere()