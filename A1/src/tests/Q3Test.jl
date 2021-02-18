using Revise

using Memento
using Plots
using Printf
using Test

using ForwardDiff


includet("../A1_module/A1Module.jl")
using .A1Module

setlevel!(getlogger("Main"), "debug")
setlevel!(getlogger(A1Module), "debug")
Memento.config!("debug")

const TIME_DATA = [0; 1; 2; 3; 4; 5; 6; 7; 8; 9]
const Y_DATA = [1.75; 1.65; 1.56; 1.49; 1.43; 1.37; 1.33; 1.29; 1.26; 1.23]

N_f_eval = 0
N_grad_f_eval = 0

function Q3Model(t, params)
    return @. params[1] + params[2]*exp(params[4]*t) + params[3]*exp(params[5]*t)
end

function Q3_Grad_Model(t, params)
    gradient = [1;
                exp(params[4]*t);
                exp(params[5]*t);
                params[2]*t*exp(params[4]*t);
                params[3]*t*exp(params[5]*t) ]
    return gradient
end


function Q3SumSquaredError(params, data_t, data_y)
    global N_f_eval += 1

    squared_errors = (data_y - Q3Model(data_t, params)).^2
    return sum(squared_errors)
end

function Q3_Grad_SumSquaredError(params, data_t, data_y)
    global N_grad_f_eval += 1
    gradient = zeros(5)

    for (t_i, y_i) in zip(data_t, data_y)
        gradient += 2*(y_i - Q3Model(t_i, params))*Q3_Grad_Model(t_i, params)
    end
    @assert size(gradient) == size(params)

    return gradient
end


objective_function = params -> Q3SumSquaredError(params, TIME_DATA, Y_DATA)
grad_objective_function = params -> Q3_Grad_SumSquaredError(params, TIME_DATA, Y_DATA)
autodiff_grad_objective_function = params -> ForwardDiff.gradient(objective_function, params)

N_f_eval = 0
N_grad_f_eval = 0
result, history = Q2SteepestDescent(objective_function, grad_objective_function, [.1, .1, .1, .1, .1], 1e-4;  linesearch_method = "SwannsBracketingMethod")

best_params = result

t_bestfit = 0:0.1:9
y_bestfit = Q3Model(t_bestfit, best_params)
plot(t_bestfit, y_bestfit, label="Best Fit")
scatter!(TIME_DATA, Y_DATA, label="Data Points")

begin
    is, points = get(history)
    errors = []
    for (i, current_params) in enumerate(history)
        error = Q3SumSquaredError(current_params, TIME_DATA, Y_DATA)
        push!(errors, error)
    end
    plot(is, errors, yscale=:log10, shape=:circle, markersize=3)
    title!("Sum Squared Error vs Gradient Descent Iteration")
    xlabel!("Number of Gradient Descent Iterations")
    ylabel!("Sum Squared Error")
end