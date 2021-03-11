include("makeplots.jl")

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

println("\nHookeJeeves")
@time evaluateHookeJeeves()

println("\nNelder Mead")
@time evaluateNelderMead()

println("\nOriginal Newtons Method")
@time evaluateOriginalNewtonsMethod()

println("\nModified Newtons Method with Levenberg Marquardt")
@time evaluateModifiedNewtonsWithLM()