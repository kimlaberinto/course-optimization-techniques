include("makeplots.jl")

println("\nHookeJeeves")
@time evaluateHookeJeeves()

println("\nGradient Descent")
@time evaluateGradientDescent()

println("\nPowell Conjugate Gradient Descent")
@time evaluatePowellConjugateGradient()

println("\nConjugate Gradient (Fletcher-Reeves)")
@time evaluateConjugateGradientFletcherReeves()

println("\nConjugate Gradient (Hestenes-Stiefel)")
@time evaluateConjugateGradientHestenesStiefel()

println("\nConjugate Gradient (Polak-Ribiere)")
@time evaluateConjugateGradientPolakRibiere()

println("\nOriginal Newtons Method")
@time evaluateOriginalNewtonsMethod()