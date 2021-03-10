module A2Module

#Importing A1Module 
include("A1Module.jl") # Module written by me, Kim Laberinto
using .A1Module: SwannsBracketingMethod, GoldenSectionSearch

# Other Imports
using LinearAlgebra #External module for taking norm
using Memento #Invenia Module for logging
using Printf #External module for formatting strings
using ValueHistories #External Package for keeping track of values

# Exports
export conjugateGradient

# Set up Memento Logger
const LOGGER = getlogger(@__MODULE__)
function __init__()
    Memento.register(LOGGER)
end

function conjugateGradient(f::Function, grad_f::Function, x_0::Array, 
        g_tol::T, k_max::Integer, n_resetsearchdir::Integer; method = "",
        tol_for_linesearch = 1e-3) where T <: Real
    info(LOGGER, "Entering Conjugate Gradient ($method)")

    k_current = 0;
    x_current = x_0;
    num_resetsearchdir = 0;

    g_old = grad_f(x_current);
    g_new = g_old;
    
    s_old = g_old;
    s_new = g_old;

    history = MVHistory()
    push!(history, :x_current, 0, x_current)
    push!(history, :grad_norm, 0, norm(g_new))
    push!(history, :num_resets, 0, num_resetsearchdir)

    while (norm(g_new) > g_tol) && (k_current < k_max)
        for i in 1:n_resetsearchdir
            k_current += 1 # NOTE: k moved to inner loop to track all iterations
            if i == 1
                debug(LOGGER, "Reseting Search Direction at k_current=$k_current")
                s_new = -1 * g_new
                num_resetsearchdir += 1
            else
                if method == "FletcherReeves"
                    gamma = norm(g_new)^2 / norm(g_old)^2
                elseif method == "HestenesStiefel"
                    error("Undefined method '$method' for Conjugate-Gradient Descent")
                elseif method == "PolakRibiere"
                    error("Undefined method '$method' for Conjugate-Gradient Descent")
                else
                    error("Undefined method '$method' for Conjugate-Gradient Descent")
                end

                s_new = -g_new + gamma*s_old;
            end

            debug(LOGGER, "Entering Line Search...")
            full_Nd_minimizer, _, _ = A1Module.Q1LineSearch(f, s_new, x_current, 
                tol_for_linesearch; linesearch_method = "SwannsBracketingMethod")
            # full_Nd_minimizer is already x_current + lambda*s_new
            x_current = full_Nd_minimizer;
            g_old = g_new
            g_new = grad_f(x_current);

            s_old = s_new;

            push!(history, :x_current, k_current, x_current)
            push!(history, :grad_norm, k_current, norm(g_new))
            push!(history, :num_resets, k_current, num_resetsearchdir)
        end
    end

    info(LOGGER, "Exiting Conjugate Gradient ($method)")
    return x_current, history
end

end
