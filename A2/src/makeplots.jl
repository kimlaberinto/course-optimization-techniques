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

function generatePlot_LossVsIterations(array_of_histories::Array{MVHistory{History}}, 
    array_of_labels::Array{String},
    symbol_to_get::Symbol)

    @assert length(array_of_histories) == length(array_of_labels)

    resultant_plot = plot()

    for (label, historyofhistories) in zip(array_of_labels, array_of_histories)
        is, xs = get(historyofhistories, symbol_to_get)
        
        errors = []
        for (i, x) in zip(is, xs)
            error = Rosenbrock5D(x)
            push!(errors, error)
        end

        plot!(resultant_plot, is, errors, label=label,
            yscale=:log10, lw=3, shape = :circle, markersize=3)
    end

    return resultant_plot
end

function makeDataDict(initial_vector, final_vector, final_loss; 
        N_f_evals = 0, N_grad_evals = 0, N_hessian_evals = 0)

    return Dict(
        "initial_vector" => initial_vector,
        "final_vector" => final_vector,
        "final_loss" => final_loss,
        "N_f_evals" => N_f_evals,
        "N_grad_evals" => N_grad_evals,
        "N_hessian_evals" => N_hessian_evals
    )

end

function onerunHookeJeeves(x_0::Array{Float64})
    initial_delta = 1.;
    final_delta = 1.e-3;
    orthogonal_directions = [[1., 0., 0., 0., 0.],
                             [0., 1., 0., 0., 0.],
                             [0., 0., 1., 0., 0.],
                             [0., 0., 0., 1., 0.],
                             [0., 0., 0., 0., 1.]];
    global N_f_evals = 0;
    best_result, history = HookeJeeves(Rosenbrock5D, x_0, initial_delta, final_delta, 
        orthogonal_directions)
    final_loss = Rosenbrock5D(best_result)

    data_dict = makeDataDict(x_0, best_result, final_loss; N_f_evals = N_f_evals)

    return data_dict, best_result, history
end

function evaluateHookeJeeves()

    array_of_inits = [[0., 0., 0., 0., 0.],
    (rand(5).*2 .- 1).*2,
    (rand(5).*2 .- 1).*2];

    array_of_labels = ["Initial Vector $i" for i in 1:length(array_of_inits)];
    array_of_trials_dicts = Array{Dict}(undef, length(array_of_inits));
    array_of_histories = Array{MVHistory{History}}(undef, length(array_of_inits));

    for (i, (label, x_0)) in enumerate(zip(array_of_labels, array_of_inits))
        data_dict, best_result, history = onerunHookeJeeves(x_0)
        @show typeof(history)
        @show typeof(array_of_histories)
        @show data_dict
        @show best_result

        array_of_trials_dicts[i] = Dict(label => data_dict)
        array_of_histories[i] = history
    end

    @show length(array_of_histories)
    generatePlot_LossVsIterations(array_of_histories, array_of_labels, :x_1)
    plot!()
end


