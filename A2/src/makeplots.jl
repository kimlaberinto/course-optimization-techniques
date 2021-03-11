include("A1Module.jl")
include("A2Module.jl")
include("objectivefunction.jl")

# Self-written modules import
using .A1Module: HookeJeeves, Q2SteepestDescent
using .A2Module
using .objectivefunctionModule: NDRosenbrock, autodiffGradientNDRosenbrock,
    autodiffHessianNDRosenbrock

# Useful external modules
using LinearAlgebra: cond
using Memento
using OrderedCollections
using Plots
using ValueHistories: MVHistory, History
import YAML

# Suppress Memento from inner modules
setlevel!(getlogger(A1Module), "not_set")
setlevel!(getlogger(A2Module.A1Module), "not_set")

# Define general settings (N = 5 dimensional rosenbrock)
const test_initial_point = zeros(5);

#Generated Initial Vectors 2 to 5 are generated using generate_random_inits.jl
const array_of_inits = [[ 0.00,  0.00,  0.00,  0.00,  0.00],
                        [ 0.36,  1.18, -1.01, -1.73, -0.90],
                        [ 1.07,  1.42,  0.32,  1.83,  0.61],
                        [ 0.26, -1.20,  0.60,  0.59, -1.77],
                        [-0.16, -0.81, -1.96, -1.55,  1.37]]; 

# Define objective functions to use
global N_f_evals = 0
global N_grad_evals = 0
global N_hessian_evals = 0

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
        N_f_evals = 0, N_grad_evals = 0, N_hessian_evals = 0, N_linsys_solves = 0)

    return OrderedDict(
        "initial_vector" => initial_vector,
        "final_vector" => final_vector,
        "final_loss" => final_loss,
        "N_f_evals" => N_f_evals,
        "N_grad_evals" => N_grad_evals,
        "N_hessian_evals" => N_hessian_evals,
        "N_linsys_solves" => N_linsys_solves
    )

end

function onerunGradientDescent(x_0::Array{Float64})
    tol = 1e-4;
    global N_f_evals = 0;
    global N_grad_evals = 0;
    best_result, history = A1Module.Q2SteepestDescent(Rosenbrock5D, 
        GradRosenbrock5D, x_0, tol;
        linesearch_method="SwannsBracketingMethod")
    final_loss = Rosenbrock5D(best_result)
    
    data_dict = makeDataDict(x_0, best_result, final_loss; 
        N_f_evals = N_f_evals,
        N_grad_evals = N_grad_evals)

    return data_dict, best_result, history
end

function onerunPowellConjugateGradient(x_0::Array{T}) where T <: Real
    tol = 1e-6
    linesearch_tol = 1e-3
    max_iter = 10000;
    global N_f_evals = 0;
    best_result, history = powellsConjugateGradientMethod(Rosenbrock5D, x_0, 
        tol; max_iter = max_iter, linesearch_tol=linesearch_tol)
    final_loss = Rosenbrock5D(best_result)

    data_dict = makeDataDict(x_0, best_result, final_loss;
        N_f_evals = N_f_evals)

    return data_dict, best_result, history
end

function onerunConjugateGradient(x_0::Array{T}, method::String) where T <: Real
    tol_for_linesearch = 1e-3;
    g_tol = 1e-4
    k_max = 10000;
    n_resetsearchdir = 10;

    global N_f_evals = 0;
    global N_grad_evals = 0;
    best_result, history = conjugateGradient(Rosenbrock5D, GradRosenbrock5D, x_0, 
        g_tol, k_max, n_resetsearchdir; 
        method=method, tol_for_linesearch=tol_for_linesearch)
    final_loss = Rosenbrock5D(best_result)
    
    data_dict = makeDataDict(x_0, best_result, final_loss; 
        N_f_evals = N_f_evals,
        N_grad_evals = N_grad_evals)

    return data_dict, best_result, history
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

function onerunOriginalNewtonsMethod(x_0::Array{T}) where T <: Real
    g_tol = 1e-3
    max_iter = 1000
    global N_f_evals = 0
    global N_grad_evals = 0
    global N_hessian_evals = 0
    best_result, history, num_linsys_solves = originalNewtonsMethod(GradRosenbrock5D, 
        HessianRosenbrock5D, x_0; g_tol = g_tol, max_iter = max_iter)
    final_loss = Rosenbrock5D(best_result)

    data_dict = makeDataDict(x_0, best_result, final_loss; 
        N_f_evals = N_f_evals,
        N_grad_evals = N_grad_evals,
        N_hessian_evals = N_hessian_evals,
        N_linsys_solves = num_linsys_solves)

    return data_dict, best_result, history
end


function evaluateGradientDescent()
    array_of_labels = ["Initial Vector $i" for i in 1:length(array_of_inits)];
    array_of_trials_dicts = Array{OrderedDict}(undef, length(array_of_inits));
    array_of_histories = Array{MVHistory{History}}(undef, length(array_of_inits));

    for (i, (label, x_0)) in enumerate(zip(array_of_labels, array_of_inits))
        data_dict, best_result, history = onerunGradientDescent(x_0)
        # @show typeof(history)
        # @show typeof(array_of_histories)
        # @show data_dict
        # @show best_result

        array_of_trials_dicts[i] = OrderedDict(label => data_dict)
        array_of_histories[i] = history
    end
    all_trial_dicts = merge(array_of_trials_dicts...)
    YAML.write_file("assets/GradientDescent_TrialOutputs.yml", all_trial_dicts)


    plot_losses = generatePlot_LossVsIterations(array_of_histories, array_of_labels, :Nd_point)
    xlabel!(plot_losses, "Number of Gradient Descent Iterations")
    ylabel!(plot_losses, "Loss")
    title!(plot_losses, "Loss vs Iterations - Gradient Descent")
    savefig(plot_losses, "assets/GradientDescentLossPlot.svg")
end

function evaluatePowellConjugateGradient()
    array_of_labels = ["Initial Vector $i" for i in 1:length(array_of_inits)];
    array_of_trials_dicts = Array{OrderedDict}(undef, length(array_of_inits));
    array_of_histories = Array{MVHistory{History}}(undef, length(array_of_inits));

    for (i, (label, x_0)) in enumerate(zip(array_of_labels, array_of_inits))
        data_dict, best_result, history = onerunPowellConjugateGradient(x_0)

        array_of_trials_dicts[i] = OrderedDict(label => data_dict)
        array_of_histories[i] = history
    end

    all_trial_dicts = merge(array_of_trials_dicts...)
    YAML.write_file("assets/PowellConjugateGradient_TrialOutputs.yml", all_trial_dicts)

    plot_losses = generatePlot_LossVsIterations(array_of_histories, array_of_labels, :x_current)
    xlabel!(plot_losses, "Number of Iterations")
    ylabel!(plot_losses, "Loss")
    title!(plot_losses, "Loss vs Iterations\nPowell Conjugate Gradient")
    savefig(plot_losses, "assets/PowellConjugateGradient_LossPlot.svg")
end

function evaluateConjugateGradientFletcherReeves()
    array_of_labels = ["Initial Vector $i" for i in 1:length(array_of_inits)];
    array_of_trials_dicts = Array{OrderedDict}(undef, length(array_of_inits));
    array_of_histories = Array{MVHistory{History}}(undef, length(array_of_inits));

    for (i, (label, x_0)) in enumerate(zip(array_of_labels, array_of_inits))
        data_dict, best_result, history = onerunConjugateGradient(x_0, "FletcherReeves")
        # @show typeof(history)
        # @show typeof(array_of_histories)
        # @show data_dict
        # @show best_result

        array_of_trials_dicts[i] = OrderedDict(label => data_dict)
        array_of_histories[i] = history
    end

    all_trial_dicts = merge(array_of_trials_dicts...)
    YAML.write_file("assets/ConjugateGradientFletcherReeves_TrialOutputs.yml", all_trial_dicts)

    plot_losses = generatePlot_LossVsIterations(array_of_histories, array_of_labels, :x_current)
    xlabel!(plot_losses, "Number of Iterations/Updates")
    ylabel!(plot_losses, "Loss")
    title!(plot_losses, "Loss vs Iterations\nConjugate Gradient (Fletcher-Reeves)")
    savefig(plot_losses, "assets/ConjugateGradientFletcherReeves_LossPlot.svg")
end

function evaluateConjugateGradientHestenesStiefel()
    array_of_labels = ["Initial Vector $i" for i in 1:length(array_of_inits)];
    array_of_trials_dicts = Array{OrderedDict}(undef, length(array_of_inits));
    array_of_histories = Array{MVHistory{History}}(undef, length(array_of_inits));

    for (i, (label, x_0)) in enumerate(zip(array_of_labels, array_of_inits))
        data_dict, best_result, history = onerunConjugateGradient(x_0, "HestenesStiefel")
        # @show typeof(history)
        # @show typeof(array_of_histories)
        # @show data_dict
        # @show best_result

        array_of_trials_dicts[i] = OrderedDict(label => data_dict)
        array_of_histories[i] = history
    end

    all_trial_dicts = merge(array_of_trials_dicts...)
    YAML.write_file("assets/ConjugateGradientHestenesStiefel_TrialOutputs.yml", all_trial_dicts)

    plot_losses = generatePlot_LossVsIterations(array_of_histories, array_of_labels, :x_current)
    xlabel!(plot_losses, "Number of Iterations/Updates")
    ylabel!(plot_losses, "Loss")
    title!(plot_losses, "Loss vs Iterations\nConjugate Gradient (Hestenes-Stiefel)")
    savefig(plot_losses, "assets/ConjugateGradientHestenesStiefel_LossPlot.svg")
end

function evaluateConjugateGradientPolakRibiere()
    array_of_labels = ["Initial Vector $i" for i in 1:length(array_of_inits)];
    array_of_trials_dicts = Array{OrderedDict}(undef, length(array_of_inits));
    array_of_histories = Array{MVHistory{History}}(undef, length(array_of_inits));

    for (i, (label, x_0)) in enumerate(zip(array_of_labels, array_of_inits))
        data_dict, best_result, history = onerunConjugateGradient(x_0, "PolakRibiere")
        # @show typeof(history)
        # @show typeof(array_of_histories)
        # @show data_dict
        # @show best_result

        array_of_trials_dicts[i] = OrderedDict(label => data_dict)
        array_of_histories[i] = history
    end

    all_trial_dicts = merge(array_of_trials_dicts...)
    YAML.write_file("assets/ConjugateGradientPolakRibiere_TrialOutputs.yml", all_trial_dicts)

    plot_losses = generatePlot_LossVsIterations(array_of_histories, array_of_labels, :x_current)
    xlabel!(plot_losses, "Number of Iterations/Updates")
    ylabel!(plot_losses, "Loss")
    title!(plot_losses, "Loss vs Iterations\nConjugate Gradient (Polak-Ribière)")
    savefig(plot_losses, "assets/ConjugateGradientPolakRibiere_LossPlot.svg")
end

function evaluateHookeJeeves()
    array_of_labels = ["Initial Vector $i" for i in 1:length(array_of_inits)];
    array_of_trials_dicts = Array{OrderedDict}(undef, length(array_of_inits));
    array_of_histories = Array{MVHistory{History}}(undef, length(array_of_inits));

    for (i, (label, x_0)) in enumerate(zip(array_of_labels, array_of_inits))
        data_dict, best_result, history = onerunHookeJeeves(x_0)
        # @show typeof(history)
        # @show typeof(array_of_histories)
        # @show data_dict
        # @show best_result

        array_of_trials_dicts[i] = OrderedDict(label => data_dict)
        array_of_histories[i] = history
    end

    all_trial_dicts = merge(array_of_trials_dicts...)
    YAML.write_file("assets/HookeJeeves_TrialOutputs.yml", all_trial_dicts)


    plot_losses = generatePlot_LossVsIterations(array_of_histories, array_of_labels, :x_1)
    xlabel!(plot_losses, "Number of Hooke-Jeeves Outer-loop Iterations")
    ylabel!(plot_losses, "Loss")
    title!(plot_losses, "Loss vs Iterations - Hooke-Jeeves")
    savefig(plot_losses, "assets/HookeJeevesLossPlot.svg")
end

function evaluateOriginalNewtonsMethod()
    array_of_labels = ["Initial Vector $i" for i in 1:length(array_of_inits)];
    array_of_trials_dicts = Array{OrderedDict}(undef, length(array_of_inits));
    array_of_histories = Array{MVHistory{History}}(undef, length(array_of_inits));

    for (i, (label, x_0)) in enumerate(zip(array_of_labels, array_of_inits))
        data_dict, best_result, history = onerunOriginalNewtonsMethod(x_0)

        array_of_trials_dicts[i] = OrderedDict(label => data_dict)
        array_of_histories[i] = history
    end

    all_trial_dicts = merge(array_of_trials_dicts...)
    YAML.write_file("assets/OriginalNewtonsMethod_TrialOutputs.yml", all_trial_dicts)

    plot_losses = generatePlot_LossVsIterations(array_of_histories, array_of_labels, :x_current)
    xlabel!(plot_losses, "Number of Iterations")
    ylabel!(plot_losses, "Loss")
    title!(plot_losses, "Loss vs Iterations\nOriginal Newtons Method")
    savefig(plot_losses, "assets/OriginalNewtonsMethod_LossPlot.svg")

    begin
        plot_condnum_matrices = plot()
        for (label, historyofhistories) in zip(array_of_labels, array_of_histories)
            is, matrices = get(historyofhistories, :hessian_current)
            condnums = []
            for (i, M) in zip(is, matrices)
                condnum = cond(M)
                push!(condnums, condnum)
            end
            plot!(plot_condnum_matrices, is, condnums, label=label,
                yscale=:log10, lw=3, shape = :circle, markersize=3, legend=:bottomright)
        end
        xlabel!(plot_condnum_matrices, "Number of Iterations")
        ylabel!(plot_condnum_matrices, "Condition Number of Hessian")
        title!(plot_condnum_matrices, "Hessian Condition Number (log scale) vs Iterations\nOriginal Newtons Method")
        savefig(plot_condnum_matrices, "assets/OriginalNewtonsMethod_ConditionNumberHessianPlot.svg")
    end
end