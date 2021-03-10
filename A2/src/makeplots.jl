include("A1Module.jl")
include("A2Module.jl")
include("objectivefunction.jl")

# Self-written modules import
using .A1Module: HookeJeeves, Q2SteepestDescent
using .A2Module
using .objectivefunctionModule: NDRosenbrock, autodiffGradientNDRosenbrock,
    autodiffHessianNDRosenbrock

# Useful external modules
using Plots
using ValueHistories: MVHistory, History

# Define general settings (N = 5 dimensional rosenbrock)
const test_initial_point = zeros(5);

# Define objective functions to use
global N_f_evals = 0
global N_grad_evals = 0
global N_hessian_evals = 0
global N_linearsystemsolves = 0 #TODO

function Rosenbrock5D(x::Array{T}) where T <: Real
    global N_f_evals += 1
    return NDRosenbrock(5, x)
end

function GradRosenbrock5D(x::Array{T}) where T <: Real
    global N_grad_evals += 1
    return autodiffGradientNDRosenbrock(5, x)
end

function HessianRosenbrock5D(x::Array{T}) where T <: Real
    global N_hessian_evals += 1
    return autodiffHessianNDRosenbrock(5, x)
end

function generatePlot_LossVsIterations(historyofhistories::MVHistory{History}, 
    symbol_to_get::Symbol)

    resultant_plot = plot()
    is, xs = get(historyofhistories, symbol_to_get)
    
    errors = []
    for (i, x) in zip(is, xs)
        error = Rosenbrock5D(x)
        push!(errors, error)
    end
    
    resultant_plot = plot(is, errors, yscale=:log10, lw=3, shape = :circle, markersize=3)
end

function evaluateHookeJeeves()
    x_0 = [0., 0., 0., 0., 0.];
    initial_delta = 1.;
    final_delta = 1.e-3;
    orthogonal_directions = [[1., 0., 0., 0., 0.],
                             [0., 1., 0., 0., 0.],
                             [0., 0., 1., 0., 0.],
                             [0., 0., 0., 1., 0.],
                             [0., 0., 0., 0., 1.]]
    best_result, history = HookeJeeves(Rosenbrock5D, x_0, initial_delta, final_delta, 
        orthogonal_directions)

    generatePlot_LossVsIterations(history, :x_1)
    plot!()

    global N_f_evals
    @show N_f_evals
end