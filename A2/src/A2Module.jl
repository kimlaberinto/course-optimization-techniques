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
export secantLineSearch
export powellsConjugateGradientMethod

# Set up Memento Logger
const LOGGER = getlogger(@__MODULE__)
function __init__()
    Memento.register(LOGGER)
end

function secantLineSearch(grad_f::Function, x_0::Array, d::Array, linesearch_tol::T;
        max_iter::Integer = 100) where T <: Real
    debug(LOGGER, "Entering Secant Line Search")
    alpha_current = 0.0;
    alpha = 0.1; #Larger initial alpha. No more NaNs
    dphi_zero = grad_f(x_0)' * d
    dphi_current = dphi_zero

    i = 0
    while abs(dphi_current) > linesearch_tol*abs(dphi_zero)
        alpha_old = alpha_current;
        alpha_current = alpha;
        dphi_old = dphi_current;
        dphi_current = grad_f(x_0 .+ alpha_current.*d)' * d;
        alpha = (dphi_current * alpha_old - dphi_old * alpha_current) / (dphi_current - dphi_old);
        i += 1

        if i >= max_iter && abs(dphi_current) < abs(dphi_zero)
            debug(LOGGER, "Secant Line Search Terminating with i=$i")
            break
        end
    end

    full_Nd_point = x_0 .+ alpha_current.*d
    debug(LOGGER, "Exiting Secant Line Search")
    return full_Nd_point
end

function powellsConjugateGradientMethod(f::Function, x_0::Array, tol::T; 
    max_iter::Integer = 1000, linesearch_tol = 1e-3) where T <: Real
    info(LOGGER, "Entering Powells Conjugate Gradient")

    search_dir_array = Array{Array}(undef, length(x_0)+1)
    for i in 1:length(x_0)
        search_dir_array[i] = zeros(length(x_0))
        search_dir_array[i][i] = 1
    end

    full_Nd_minimizer, _, _ = A1Module.Q1LineSearch(f, search_dir_array[length(x_0)], x_0, linesearch_tol;
        linesearch_method = "SwannsBracketingMethod")
    X = full_Nd_minimizer;
    C = false;
    k = 0;

    history = MVHistory()
    push!(history, :x_current, 0, x_0)

    while C == false
        Y = X;
        k += 1;
        for i in 1:length(x_0)
            #Keep updating X using line searches in the s_i directions
            # @show search_dir_array[i]
            # @show X
            full_Nd_minimizer, _, _ = A1Module.Q1LineSearch(f, search_dir_array[i], X, 
                linesearch_tol; linesearch_method = "SwannsBracketingMethod")
            X = full_Nd_minimizer;
        end
        search_dir_array[end] = X .- Y;

        full_Nd_minimizer, _, _ = A1Module.Q1LineSearch(f, search_dir_array[end], X, 
            linesearch_tol; linesearch_method = "SwannsBracketingMethod")
        X = full_Nd_minimizer;

        f_X = f(X)
        f_Y = f(Y)
        if k > max_iter || (abs(f_X - f_Y) / max(abs(f_X), 1e-10)) < tol
            C = true
        else
            for i in 1:length(x_0)
                search_dir_array[i] = search_dir_array[i+1]
            end
        end

        history = MVHistory()
        push!(history, :x_current, k, X)
    end

    info(LOGGER, "Exiting Powells Conjugate Gradient")
    return X, history
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
                    gamma = (g_old' * g_new) / (g_old' * s_old)
                elseif method == "PolakRibiere"
                    gamma = (g_old' * g_new) / norm(g_old)^2
                else
                    error("Undefined method '$method' for Conjugate-Gradient Descent")
                end

                s_new = -g_new + gamma*s_old;
            end

            debug(LOGGER, "Entering Line Search...")
            # full_Nd_minimizer, _, _ = A1Module.Q1LineSearch(f, s_new, x_current, 
            #    tol_for_linesearch; linesearch_method = "SwannsBracketingMethod")
            full_Nd_minimizer = secantLineSearch(grad_f, x_current, s_new, tol_for_linesearch)
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
