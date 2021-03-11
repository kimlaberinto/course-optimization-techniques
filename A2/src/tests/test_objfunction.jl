using Test

include("../objectivefunction.jl")
using .objectivefunctionModule

@testset "Acceptable Inputs of Rosenbrock" begin
    @test_throws AssertionError NDRosenbrock(0, ones(2))
    @test_throws AssertionError NDRosenbrock(1, ones(2))
end

@testset "Zeros of Rosenbrock" begin
    for N in 2:10
        @test NDRosenbrock(N, ones(N)) == 0
    end
end

@testset "Gradients" begin
    @testset "Gradients at origin" begin
        for N in 2:10
            @test isapprox(
                analyticGradientNDRosenbrock(N, zeros(N)),
                autodiffGradientNDRosenbrock(N, zeros(N)))
        end
    end
end

# @testset "Hessians" begin
#     @testset "Hessians at origin" begin
#         for N in 2:10
#             @test isapprox(
#                 analyticHessianNDRosenbrock(N, zeros(N)),
#                 autodiffHessianNDRosenbrock(N, zeros(N)))
#         end
#     end
# end