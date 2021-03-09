using Test

include("../objectivefunction.jl")

@testset "Acceptable Inputs of Rosenbrock" begin
    @test_throws AssertionError NDRosenbrock(0, ones(2))
    @test_throws AssertionError NDRosenbrock(1, ones(2))
end

@testset "Zeros of Rosenbrock" begin
    for i in 2:10
        @test NDRosenbrock(i, ones(i)) == 0
    end
end