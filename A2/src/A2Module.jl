module A2Module

#Importing A1Module 
include("A1Module.jl") # Module written by me, Kim Laberinto
using .A1Module: SwannsBracketingMethod, GoldenSectionSearch

# Other Imports
using LinearAlgebra #External module for taking norm, condition number, and identity
using Memento #Invenia Module for logging
using Printf #External module for formatting strings
using ValueHistories #External Package for keeping track of values

# Exports
export conjugateGradient
export secantLineSearch
export powellsConjugateGradientMethod
export nelderMeadSimplexSearch
export originalNewtonsMethod
export modifiedNewtonsWithLMMethod

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
            try 
                full_Nd_minimizer, _, _ = A1Module.Q1LineSearch(f, search_dir_array[i], X, 
                    linesearch_tol; linesearch_method = "SwannsBracketingMethod")
                X = full_Nd_minimizer;
            catch e
                warn(LOGGER, "Caught an error during Powells Conjugate Gradient with x_0=$x_0. Aborting and returning history. $e")
                push!(history, :x_current, k, X)
                return X, history
            end
        end
        search_dir_array[end] = X .- Y;

        try
            full_Nd_minimizer, _, _ = A1Module.Q1LineSearch(f, search_dir_array[end], X, 
                linesearch_tol; linesearch_method = "SwannsBracketingMethod")
            X = full_Nd_minimizer;
        catch e
            warn(LOGGER, "Caught an error during Powells Conjugate Gradient with x_0=$x_0. Aborting and returning history. $e")
            push!(history, :x_current, k, X)
            return X, history
        end

        f_X = f(X)
        f_Y = f(Y)
        if k > max_iter || (abs(f_X - f_Y) / max(abs(f_X), 1e-10)) < tol
            C = true
        else
            for i in 1:length(x_0)
                search_dir_array[i] = search_dir_array[i+1]
            end
        end

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

function originalNewtonsMethod(grad_f::Function, hessian_f::Function, 
        x_0::Array{T}; g_tol::T = 1e-3, max_iter::Integer = 1000) where T <: Real
    info(LOGGER, "Entering Original Newtons Method")

    k_current = 0;
    num_linsys_solves = 0;

    x_current = x_0;

    g_current = grad_f(x_current);

    history = MVHistory()
    push!(history, :x_current, 0, x_current)
    push!(history, :g_current, 0, g_current)

    while norm(g_current) > g_tol && k_current < max_iter
        k_current += 1

        num_linsys_solves += 1
        hessian_current = hessian_f(x_current)
        g_current = grad_f(x_current)
        d = hessian_current \ (-g_current)
        x_current += d

        push!(history, :x_current, k_current, x_current)
        push!(history, :g_current, k_current, g_current)
        push!(history, :hessian_current, k_current, hessian_current)
    end

    info(LOGGER, "Exiting Original Newtons Method")
    return x_current, history, num_linsys_solves
end

function modifiedNewtonsWithLMMethod(grad_f::Function, hessian_f::Function,
    x_0::Array{T}; linesearch_tol::T = 1e-3, mu_param::T = 1.0, g_tol::T = 1e-3, 
    max_iter::Integer = 1000) where T <: Real
    info(LOGGER, "Entering Modified Newtons Method with LM")

    k_current = 0;
    num_linsys_solves = 0;

    x_current = x_0;

    g_current = grad_f(x_current);

    history = MVHistory()
    push!(history, :x_current, 0, x_current)
    push!(history, :g_current, 0, g_current)

    while norm(g_current) > g_tol && k_current < max_iter
        k_current += 1

        num_linsys_solves += 1
        hessian_current = hessian_f(x_current)
        LM_matrix = (hessian_current + I*mu_param)
        g_current = grad_f(x_current)
        d = LM_matrix \ (-g_current)
        
        #Do a line search to do the update
        x_current = secantLineSearch(grad_f, x_current, d, linesearch_tol)

        push!(history, :x_current, k_current, x_current)
        push!(history, :g_current, k_current, g_current)
        push!(history, :hessian_current, k_current, hessian_current)
        push!(history, :LM_matrix, k_current, LM_matrix)
    end

    info(LOGGER, "Exiting Modified Newtons Method with LM")
    return x_current, history, num_linsys_solves

end

# https://c.mql5.com/31/43/garch-improved-nelder-mead-mt4-screen-9584.png
function nelderMeadSimplexSearch(f::Function, x_0::Array{T}, 
    initial_sidelength::T; max_iter::Integer = 1000, stuck_max::Integer = 10,
    stuck_coef::T = 0.5) where T <: AbstractFloat
    info(LOGGER, "Entering Nelder Mead")
    current_vertices = generateSimplex(x_0, initial_sidelength)
    @assert length(current_vertices) == (length(x_0) + 1)

    k_current = 0;

    x_l_old = x_0;
    stuck_counter = 0;

    current_sidelength = initial_sidelength

    history = MVHistory()
    push!(history, :x_best, 0, x_0)

    ALPHA = 1;
    BETA = 0.5;
    GAMMA = 2;
    while k_current < max_iter
        k_current += 1
        #Evaluate and sort into ascending order.
        #current_vertices[1] will be the best (lowest) vertex
        f_evals = map(f, current_vertices)
        indices_for_sorting = sortperm(f_evals)
        current_vertices = current_vertices[indices_for_sorting]
        f_evals = f_evals[indices_for_sorting]

        x_h = current_vertices[end] #Highest (to replace)
        f_h = f_evals[end]

        x_g = current_vertices[end-1] #Second Highest
        f_g = f_evals[end-1]

        x_l = current_vertices[1] #Lower
        f_l = f_evals[1]

        push!(history, :x_best, k_current, current_vertices[1])

        if x_l_old == x_l
            stuck_counter += 1
        end

        #Define centroid based on other vertices
        x_c = zeros(length(x_0))
        for (i, vertex) in enumerate(current_vertices)
            if vertex != x_h
                x_c += vertex
            end
        end
        x_c /= length(x_0)

        # Do a normal reflection
        x_r = 2*x_c - x_h
        f_r = f(x_r)

        if stuck_counter < stuck_max
            if f_l < f_r < f_g
                theta = ALPHA;
            elseif f_r < f_l
                theta = GAMMA;
            elseif f_r > f_h
                theta = -1 * BETA;
            else # f_g < f_r < f_h
                theta = BETA;
            end
        else
            info(LOGGER, "Shrinking in Nelder Mead")
            stuck_counter = 0 #reset
            for (index, x_old) in enumerate(current_vertices)
                if index >= 2 #only modify the non-x_l vertices
                    current_vertices[index] = (x_old - x_l)*stuck_coef + x_l
                end
            end

            x_l_old = x_l;
            continue #skip to next iteration
        end

        x_l_old = x_l;
        x_new = x_h + (1+ theta)*(x_c - x_h)

        # Replace the largest with x_new
        current_vertices[end] = x_new
    end

    info(LOGGER, "Exiting in Nelder Mead")
    return current_vertices[1], history
end

function generateSimplex(basePoint::Array{T}, side_length::T) where T <: AbstractFloat
    n = length(basePoint)

    a = side_length * ( (sqrt(n+1) + (n - 1)) / (n * sqrt(2)))
    b = side_length * ( (sqrt(n+1) - 1) / (n * sqrt(2)) )

    list_of_points = Array{Array{T}}(undef, n+1)
    for i in 1:n #Generate each point
        newPoint = Array{T}(undef, n)
        for j in 1:n # Define coordinates of new point
            if i == j
                newPoint[j] = basePoint[j] + a;
            else
                newPoint[j] = basePoint[j] + b;
            end
        end
        list_of_points[i] = newPoint
    end
    list_of_points[end] = basePoint

    @assert length(list_of_points) == (n + 1)
    return list_of_points
end

end
