using Test
using Memento

include("../A2Module.jl")
using .A2Module

include("../objectivefunction.jl")
using .objectivefunctionModule: NDRosenbrock, autodiffGradientNDRosenbrock,
    autodiffHessianNDRosenbrock

# Suppress Memento from inner modules
setlevel!(getlogger(A2Module), "not_set")

function _Rosenbrock5D(x::Array{T}) where T <: Real
    return NDRosenbrock(5, x)
end

function _GradRosenbrock5D(x::Array{T}) where T <: Real
    return autodiffGradientNDRosenbrock(5, x)
end

function _HessianRosenbrock5D(x::Array{T}) where T <: Real
    return autodiffHessianNDRosenbrock(5, x)
end

@testset "Converge to Approx. True Min from Origin" begin
    @test isapprox(conjugateGradient(_Rosenbrock5D, _GradRosenbrock5D, [0.0, 0.0, 0.0, 0.0, 0.0], 
        1e-4, 10000, 10; method="FletcherReeves", tol_for_linesearch=1e-3), 
        [1., 1., 1., 1., 1.]; atol=1e-3)
end


@testset "Throws Error with Invalid Method" begin
    @test_throws ErrorException conjugateGradient(_Rosenbrock5D, _GradRosenbrock5D, [0.0, 0.0, 0.0, 0.0, 0.0], 
        1e-4, 10000, 10; 
        method="", tol_for_linesearch=1e-3)
end