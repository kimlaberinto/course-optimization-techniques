
println("Running SwannsTests.jl...")
t = @elapsed include("SwannsTests.jl")
println("done (took $t seconds). \n")

println("Running PowellsTests.jl...")
t = @elapsed include("PowellsTests.jl")
println("done (took $t seconds). \n")

println("Running GoldenTests.jl...")
t = @elapsed include("GoldenTests.jl")
println("done (took $t seconds). \n")

println("Running Q1Test.jl...")
t = @elapsed include("Q1Test.jl")
println("done (took $t seconds). \n")

println("Running Q2Test.jl...")
t = @elapsed include("Q1Test.jl")
println("done (took $t seconds). \n")