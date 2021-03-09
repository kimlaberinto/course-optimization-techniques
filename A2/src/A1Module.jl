module A1Module

using LinearAlgebra #External module for taking norm
using Memento #Invenia Module for logging
using Printf #External module for formatting strings
using Plots
using ValueHistories #External Package for keeping track of values

export SwannsBracketingMethod
export PowellsBracketingMethod
export GoldenSectionSearch
export Q1LineSearch
export Q2SteepestDescent
export HookeJeeves

# Set up Memento Logger
const LOGGER = getlogger(@__MODULE__)
function __init__()
    Memento.register(LOGGER)
end

function SwannsBracketingMethod(f, x_0, initial_step_length; magnification = 2)
    info(LOGGER, "Inside Swann's Bracketing Method")

    #Define test points
    x_l = x_0 - abs(initial_step_length) #lower
    x_u = x_0 + abs(initial_step_length) #upper
    x_m = x_0 #middle
    

    #Do function evaluations
    debug(LOGGER, "Performing initial function evaliations")
    f_l = f(x_l) #lower
    f_u = f(x_u) #upper
    f_m = f(x_0) #middle 

    history = History(Tuple{Float64, Float64, Float64})
    push!(history, 0, convert(Tuple{Float64, Float64, Float64}, (x_l, x_u, x_m)))

    # 4 Cases.
    # 1) Keep moving right
    # 2) Keep moving left
    # 3) Initial interval is bracket
    # 4) Error (non-unimodal)

    debug(LOGGER, @sprintf "x_l x_m x_u = [%3.2f, %3.2f %3.2f]" x_l x_m x_u)
    if f_l >= f_m >= f_u
        debug(LOGGER, "Case 1: Keep Moving Right")
        i = 1
        while f_u < f_m
            debug(LOGGER, @sprintf "x_l to x_m = [%3.2f, %3.2f]" x_l x_m)

            (x_l, x_m, f_l, f_m) = (x_m, x_u, f_m, f_u)
            x_u = x_u + (magnification^i) * abs(initial_step_length)
            f_u = f(x_u)

            push!(history, i, convert(Tuple{Float64, Float64, Float64}, (x_l, x_u, x_m)))
            i += 1
        end
        # f_u is now higher than f_m
        # This now forms a valid bracketing interval

        info(LOGGER, @sprintf "Outputting x_l x_u = [%3.2f %3.2f]" x_l x_u)
        return (x_l, x_u), history

    elseif f_l <= f_m <= f_u
        debug(LOGGER, "Case 2: Keep Moving Left")
        i = 1
        while f_l < f_m
            debug(LOGGER, @sprintf "x_u to x_m = [%3.2f, %3.2f]" x_l x_m)

            (x_u, x_m, f_u, f_m) = (x_m, x_l, f_m, f_l)
            x_l = x_l - (magnification^i) * abs(initial_step_length)
            f_l = f(x_l)

            push!(history, i, convert(Tuple{Float64, Float64, Float64}, (x_l, x_u, x_m)))
            i += 1
        end
        # f_l is now higher than f_m
        # This forms a valid bracketing interval
        info(LOGGER, @sprintf "Outputting x_l x_u = [%3.2f %3.2f]" x_l x_u)
        return (x_l, x_u), history

    elseif f_l >= f_m <= f_u
        debug(LOGGER, "Case 3: Initial interval is bracket")
        info(LOGGER, @sprintf "Outputting x_l x_u = [%3.2f %3.2f]" x_l x_u)
        return (x_l, x_u), history

    else # 4) Error  f_l <= f_m >= f_u
        error(LOGGER, "Case 4: Error (non-unimodal)")
    end

end

function PowellsBracketingMethod(f, a_1, delta, delta_max)
    #TODO: Gets an error when the function is perfectly quadratic e.g x -> x^2
    info(LOGGER, "Inside Powell's Bracketing Method")

    c_1 = a_1 + delta
    F_a = f(a_1)
    F_c = f(c_1)

    if F_a > F_c
        b_1 = a_1 + 2*delta
        F_b = f(b_1)
        forward = true
    else
        (b_1, c_1, a_1) = (c_1, a_1, a_1 - delta)
        (F_b, F_c, F_a) = (F_c, F_a, f(a_1))
        forward = false
    end

    debug(LOGGER, @sprintf "Direction is forward: %s" forward)
    debug(LOGGER, @sprintf " a_1, b_1, c_1 = %s" (a_1, b_1, c_1))
    debug(LOGGER, @sprintf " F_a, F_b, F_c = %s" (F_a, F_b, F_c))

    a_current = a_1
    b_current = b_1
    c_current = c_1

    a_next = a_1
    b_next = b_1
    c_next = c_1

    history = History(Tuple{Float64, Float64, Float64})
    push!(history, 0, convert(Tuple{Float64, Float64, Float64}, (a_current, b_current, c_current)))

    N_iterations = 0
    while !(F_c < F_a && F_c < F_b)
        debug(LOGGER, @sprintf "New loop iteration...")
        debug(LOGGER, @sprintf "F_a F_b F_c = %5.3f %5.3f %5.3f" F_a F_b F_c)
        a_current = a_next
        b_current = b_next
        c_current = c_next
        N_iterations += 1

        debug(LOGGER, @sprintf "a_current, b_current, c_current = [%5.3f, %5.3f, %5.3f]" a_current b_current c_current)


        p_numerator = (c_current - b_current)*F_a + (a_current - c_current)*F_b + (b_current-a_current)*F_c
        p_denominator = ((b_current - c_current)*(c_current - a_current)*(a_current - b_current))
        p = p_numerator / p_denominator
        debug(LOGGER, @sprintf "p = %5.3f" p)

        # Cases Summarized
        # 1) Moving Forward
        # 1.1) Quadratic has maximum (p <= 0)
        # 1.2) Quadratic has no maximum
        # 1.2.1) Quadratic minimum is too far away
        # 1.2.2) Quadratic minimum is within reach
        # 2) Moving Backward
        # 2.1) Quadratic has maximum (p <= 0)
        # 2.2) Quadratic has no maximum
        # 2.2.1) Quadratic minimum is too far away
        # 2.2.2) Quadratic minimum is within reach

        if forward #Moving Forward
            debug(LOGGER, @sprintf "Moving Forward...")
            if p <= 0 #Quadratic has a maximum, move as far as possible
                debug(LOGGER, @sprintf "Quadratic has a maximum. Move as far forward as possible.")
                a_next = a_current
                b_next = b_current + delta_max
                c_next = c_current

                F_b = f(b_next)
            else
                #p > 0 concave up
                x_star_num = (b_current^2 - c_current^2)*F_a + (c_current^2 - a_current^2)*F_b + (a_current^2 - b_current^2)*F_c
                x_star_denom = (b_current - c_current)*F_a + (c_current - a_current)*F_b + (a_current - b_current)*F_c
                x_star = (1/2) * x_star_num / x_star_denom 
    
                debug(LOGGER, @sprintf "p > 0, with x_star = %5.3f" x_star)

                if (x_star - b_current) > delta_max #Quadratic minimum is too far
                    debug(LOGGER, @sprintf "Quadratic minimum is too far. Moving as far forward as possible." )
                    b_next = b_current + delta_max
                else
                    debug(LOGGER, @sprintf "Quadratic minimum is reachable. (forward)" )
                    b_next = x_star
                end

                a_next = c_current
                c_next = b_current
                
                (F_a, F_c, F_b) = (F_c, F_b, f(b_next))
            end
        else #Moving Backward
            debug(LOGGER, @sprintf "Moving Backward...")
            if p <= 0 #Quadratic has a maximum, move as far as possible
                debug(LOGGER, @sprintf "Quadratic has a maximum. Move as far back as possible.")
                a_next = b_current - delta_max
                b_next = b_current
                c_next = c_current

                F_a = f(a_next)
            else
                #p > 0 concave up
                x_star_num = (b_current^2 - c_current^2)*F_a + (c_current^2 - a_current^2)*F_b + (a_current^2 - b_current^2)*F_c
                x_star_denom = (b_current - c_current)*F_a + (c_current - a_current)*F_b + (a_current - b_current)*F_c
                x_star = (1/2) * x_star_num / x_star_denom 
    
                debug(LOGGER, @sprintf "p > 0, with x_star = %5.3f" x_star)

                if (a_current - x_star) > delta_max #Quadratic minimum is too far
                    debug(LOGGER, @sprintf "Quadratic minimum is too far. Moving as far back as possible." )
                    a_next = a_current - delta_max
                else
                    debug(LOGGER, @sprintf "Quadratic minimum is reachable. (backwards)" )
                    a_next = x_star
                end

                b_next = c_current
                c_next = a_current
                (F_b, F_c, F_a) = (F_c, F_a, f(a_next))
            end
        end

        push!(history, N_iterations, convert(Tuple{Float64, Float64, Float64}, (a_next, b_next, c_next)))

        debug(LOGGER, @sprintf "a_next b_next c_next = %5.3f %5.3f %5.3f" a_next b_next c_next)
        debug(LOGGER, @sprintf "End of loop iteration...")
    end

    a_current = a_next
    b_current = b_next
    c_current = c_next

    debug(LOGGER, @sprintf "Exiting with a_current, b_current c_current = %5.3f %5.3f %5.3f" a_current b_current c_current)
    debug(LOGGER, @sprintf "Exiting with F_a F_b F_c = %5.3f %5.3f %5.3f" F_a F_b F_c)

    return (a_current, b_current), history
end

function GoldenSectionSearch(f, a, b, tolerance)
    info(LOGGER, "Inside Golden Section Search")

    TAU = 0.618_033_988_7

    F_a = f(a)
    F_b = f(b)

    c = a + (1 - TAU)*(b - a)
    F_c = f(c)

    d = b - (1 - TAU)*(b - a)
    F_d = f(d)

    a_current = a
    b_current = b
    c_current = c
    d_current = d

    a_next = a
    b_next = b
    c_next = c
    d_next = d

    history = History(Tuple{Float64, Float64, Float64, Float64})
    push!(history, 0, convert(Tuple{Float64, Float64, Float64, Float64}, (a_current, b_current, c_current, d_current)))

    N_iterations = 0
    interval_size = abs(b - a)
    while !(interval_size < tolerance)
        a_current = a_next
        b_current = b_next
        c_current = c_next
        d_current = d_next
        N_iterations += 1

        debug(LOGGER, @sprintf "Start of loop iteration... (a_current, b_current, c_current, d_current) = (%5.3f %5.3f %5.3f %5.3f)" a_current b_current c_current d_current)
        debug(LOGGER, @sprintf "Start of loop iteration... (F_a, F_b, F_c, F_d) = (%5.3f %5.3f %5.3f %5.3f)" F_a F_b F_c F_d)

        if F_c < F_d
            debug(LOGGER, "F_c < F_d case")
            a_next = a_current
            b_next = d_current
            c_next = a_next + (1 - TAU)*(b_next - a_next)
            d_next = c_current
            
            (F_a, F_b) = (F_a, F_d)
            (F_c, F_d) = (f(c_next), F_c)

        else
            debug(LOGGER, "F_c > F_d case")
            a_next = c_current
            b_next = b_current
            c_next = d_current
            d_next = b_next - (1 - TAU)*(b_next - a_next)

            (F_a, F_b) = (F_c, F_b)
            (F_c, F_d) = (F_d, f(d_next))

        end

        push!(history, N_iterations, convert(Tuple{Float64, Float64, Float64, Float64}, (a_next, b_next, c_next, d_next)))
        interval_size = abs(b_next - a_next)
        debug(LOGGER, @sprintf "End of loop iteration... interval_size = %5.3f" interval_size)
    end

    a_current = a_next
    b_current = b_next
    c_current = c_next
    d_current = d_next

    info(LOGGER, @sprintf "Exiting Golden Section Search with a_current, b_current = %5.3f %5.3f" a_current b_current)
    info(LOGGER, @sprintf "Exiting Golden Section Search with F_a, F_b = %5.3f %5.3f" F_a F_b)
    return (a_current, b_current), history

end

function Q1LineSearch(f, d, x_0, desired_interval_size; linesearch_method = "")
    info(LOGGER, "Entering Q1LineSearch...")
    info(LOGGER, @sprintf "Entering with d = %s" d)
    info(LOGGER, @sprintf "Entering with x_0 = %s" x_0)

    one_dimensional_function = alpha -> f((x_0 .+ alpha .* d))

    bracketing_history = undef
    golden_history = undef
    a_l_smaller, a_u_smaller = undef, undef
    if linesearch_method == "SwannsBracketingMethod"
        swanns_step_length = 1 #HARDCODED
        alpha_init = 0 #HARDCODED
        (alpha_lower, alpha_upper), swanns_history = SwannsBracketingMethod(one_dimensional_function, alpha_init, swanns_step_length)
        (a_l_smaller, a_u_smaller), golden_history = GoldenSectionSearch(one_dimensional_function, alpha_lower, alpha_upper, desired_interval_size)
        bracketing_history = swanns_history
    elseif linesearch_method == "PowellsBracketingMethod"
        alpha_init = 0 #HARDCODED
        powells_delta = 1 #HARDCODED
        powells_delta_max = 16 #HARDCODED
        (alpha_lower, alpha_upper), powell_history = PowellsBracketingMethod(one_dimensional_function, alpha_init, powells_delta, powells_delta_max)
        (a_l_smaller, a_u_smaller), golden_history = GoldenSectionSearch(one_dimensional_function, alpha_lower, alpha_upper, desired_interval_size)
        bracketing_history = powell_history
    else
        error(LOGGER, "Line Search Method not recognized: $linesearch_method")
    end

    a_mid = (a_l_smaller + a_u_smaller) / 2 #along line
    debug(LOGGER, @sprintf "a_middle %5.3f" a_mid)

    full_middle_point = @. x_0 + a_mid * d #In full N-D space
    info(LOGGER, @sprintf "Exiting Q1LineSearch with middle point in N-D as %s" full_middle_point)
    return full_middle_point, bracketing_history, golden_history
end

function Q2SteepestDescent(f, grad_f, x_0, tolerance_for_1D_search; linesearch_method = "")
    info(LOGGER, "Entering Q2SteepestDescent...")
    info(LOGGER, @sprintf "Entering with x_0 = %s" x_0)

    x_0 = convert(Array{Float64, 1}, x_0)

    current_point = x_0
    next_point = x_0
    steepest_descent_direction = -1 * grad_f(x_0)

    Q2_history = MVHistory()
    push!(Q2_history, :Nd_point, 0, current_point)

    N_iterations = 0
    while !(norm(steepest_descent_direction) < 10^(-4))
        info(LOGGER, @sprintf "Q2 Start of Loop Iteration... current_point = %s" current_point)
        info(LOGGER, @sprintf "Q2 Start of Loop Iteration... steepest = %s" steepest_descent_direction)
        N_iterations += 1

        next_point, bracketing_history, golden_history = Q1LineSearch(f, steepest_descent_direction, current_point, tolerance_for_1D_search; linesearch_method = linesearch_method)
        
        steepest_descent_direction = -1 * grad_f(next_point)
        current_point = next_point

        push!(Q2_history, :Nd_point, N_iterations, current_point)
        push!(Q2_history, :bracketing_history, N_iterations, bracketing_history)
        push!(Q2_history, :golden_history, N_iterations, golden_history)
        debug(LOGGER, @sprintf "End of Loop Iteration... norm of grad %s" norm(steepest_descent_direction))
    end

    info(LOGGER, @sprintf "Exiting with next_point = %s" next_point)
    return next_point, Q2_history

end

function HookeJeeves(f, x_0, big_delta, small_delta, orthogonal_directions)
    info(LOGGER, "Entering Hooke-Jeeves...")
    @assert length(x_0) == length(orthogonal_directions)

    N = length(orthogonal_directions)
    x_i_array = zeros(length(x_0), length(orthogonal_directions)+1)

    x_0_current = x_0
    x_i_array[:, 1] = x_0
    current_big_delta = big_delta

    history = MVHistory()
    push!(history, :x_1, 0, x_i_array[:, 1])
    push!(history, :x_0_current, 0, x_0_current)
    push!(history, :current_big_delta, 0, current_big_delta)
    
    iteration_number = 0
    while !(current_big_delta < small_delta)
        iteration_number += 1
        debug(LOGGER, "Start of Iteration Loop... ($iteration_number)")


        #Exploratory Moves
        debug(LOGGER, "Performing Exploratory Moves... (delta = $current_big_delta)")
        prev_x = []
        for (index, direction) in enumerate(orthogonal_directions)
            debug(LOGGER, "Trying out direction ($index): $direction")
            debug(LOGGER, @sprintf "Trying from %s" x_i_array[:, index])
            x_current = x_i_array[:, index] .+ current_big_delta * direction

            F_x_i =  f(x_i_array[:, index])
            F_x_current_forward = f(x_current)
            if F_x_current_forward < F_x_i
                x_i_array[:, index+1] = x_current
                debug(LOGGER, "Forward is good to $x_current")
            else
                x_current = x_i_array[:, index] .- current_big_delta * direction
                F_x_current_backward = f(x_current)
                
                if F_x_current_backward < F_x_i
                    x_i_array[:, index+1] = x_current
                    debug(LOGGER, "Backward is good to $x_current")

                else
                    x_i_array[:, index+1] = x_i_array[:, index]
                    debug(LOGGER, @sprintf "Stay here at %s" x_i_array[:, index])
                    @assert (F_x_current_forward >= F_x_i)
                    @assert (F_x_current_backward >= F_x_i)


                end
            end
        end

        # Pattern Move
        debug(LOGGER, @sprintf "Performing Pattern Move Branching... with x_i_array[:, N+1] = %s" x_i_array[:, N+1])
        if f(x_i_array[:, N+1]) < f(x_0_current)
            debug(LOGGER, "Performing Pattern Move...")
            x_i_array[:, 1] = x_i_array[:, N+1] .+ (x_i_array[:, N+1] .- x_0_current)
            x_0_current = x_i_array[:, N+1]
            debug(LOGGER, @sprintf "x_i_array[:, 1] updated to %s" x_i_array[:, 1])
            debug(LOGGER, @sprintf "x_0_current updated to %s" x_0_current)
        
        elseif x_i_array[:, 1] == x_0_current
            debug(LOGGER, "Pattern move not better x_i_array[:, 1] == x_0_current")
            current_big_delta = current_big_delta / 10
            debug(LOGGER, "Reduced big_delta to $current_big_delta")
        
        else
            debug(LOGGER, "Pattern Move not better: Set x_i_array[:, 1] to x_0_current")
            x_i_array[:, 1] = x_0_current
            debug(LOGGER, @sprintf "x_i_array[:, 1] updated to %s" x_i_array[:, 1])

        end

        push!(history, :x_1, iteration_number, x_i_array[:, 1])
        push!(history, :x_0_current, iteration_number, x_0_current)
        push!(history, :x_end, iteration_number, x_i_array[:, N+1])
        push!(history, :current_big_delta, iteration_number, current_big_delta)
        debug(LOGGER, "End of Iteration Loop... ($iteration_number)")
    end

    info(LOGGER, "Exiting Hooke-Jeeves...")
    return x_i_array[:, 1], history
end

end

