using Revise

using ForwardDiff # Library for Extra Plot - Calculating Gradients and Hessians Automatically from Code (not symbolic!)
using Plots # Library for plotting
using Printf # Library for formatting strings
using LinearAlgebra # Library for calculating norm and linear algebra
using LaTeXStrings # Library for supporting LaTeX symbols in strings

using Memento #For debugging

includet("A1_module/A1Module.jl")
using .A1Module # Custom Library written for Assignment. Written by Kim Laberinto

#setlevel!(getlogger("Main"), "debug")
#setlevel!(getlogger(A1Module), "debug")
#Memento.config!("debug")

N_f_eval = 0
N_grad_f_eval = 0

## 
## SET UP
##
const TIME_DATA = [0; 1; 2; 3; 4; 5; 6; 7; 8; 9]
const Y_DATA = [1.75; 1.65; 1.56; 1.49; 1.43; 1.37; 1.33; 1.29; 1.26; 1.23]

N_f_eval = 0
N_grad_f_eval = 0

function Q3Model(t, params)
    @assert length(params) == 5
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

orthogonal_directions = [[1., 0., 0., 0., 0.],
[0., 1., 0., 0., 0.],
[0., 0., 1., 0., 0.],
[0., 0., 0., 1., 0.],
[0., 0., 0., 0., 1.]]
#orthogonal_directions = reverse(orthogonal_directions)

objective_function = params -> Q3SumSquaredError(params, TIME_DATA, Y_DATA)
const STARTING_DELTA = 0.1
const TOLERANCE_FOR_FINAL_DELTA = 1.0e-6
const STARTING_X_VECTOR = [0., 0., 0., 0., 0.]
result, history = HookeJeeves(objective_function, STARTING_X_VECTOR, STARTING_DELTA, TOLERANCE_FOR_FINAL_DELTA, orthogonal_directions)
final_param_vector = result


# Plot curve
begin
t_bestfit = 0:0.02:9
y_bestfit = Q3Model(t_bestfit, final_param_vector)
plot_curve = plot(t_bestfit, y_bestfit, label="Best Fit using Guess", legend=:topright)
end

title!("Best Fits from Optimization\nHooke-Jeeves")
xlabel!("t")
ylabel!("y")
scatter!(TIME_DATA, Y_DATA, label="Data Points")

# Plot loss curve
begin
    is, _ = get(history, :x_1)
    errors = []
    for (i, current_params) in enumerate(history, :x_1)
        error = Q3SumSquaredError(current_params, TIME_DATA, Y_DATA)
        #println("loss: $error | params: $current_params")
        push!(errors, error)
    end
    plot_loss = plot(is, errors, yscale=:log10, lw=3, label="", shape = :circle, markersize=3)
    title!("Loss vs Iteration")
end

begin
    is, current_big_deltas = get(history, :current_big_delta)
    plot_deltas = plot(is, current_big_deltas,  yscale=:log10, lw = 3, shape = :circle, markersize=3)
    title!("Current Delta vs Iteration")
end

plot_layout = @layout [a; b; c]
plot(plot_curve, plot_loss, plot_deltas, layout = plot_layout, size=(500, 700))

savefig("assets/2021-02-27_HookeJeeves_FixBug_Plot.png")
plot!() #Show plot in VSCode

println("big_delta - Starting delta: $STARTING_DELTA")
println("small_delta - Tolerance for final delta: $TOLERANCE_FOR_FINAL_DELTA")
println("Starting x vector: $STARTING_X_VECTOR")
println("Final result vector: $final_param_vector")
println(@sprintf "Final Loss: %f " Q3SumSquaredError(final_param_vector, TIME_DATA, Y_DATA))
println(@sprintf "Number of Hooke-Jeeves outerloop iterations: %s" is[end])
println(@sprintf "Number of Objective Function evaluations: %s" N_f_eval)

plot!() #Show plot in VSCode

