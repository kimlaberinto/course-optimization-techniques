module A1Module

using Memento #Invenia Module for logging
using Printf

export SwannsBracketingMethod
export PowellsBracketingMethod

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

            i += 1
        end
        # f_u is now higher than f_m
        # This now forms a valid bracketing interval

        info(LOGGER, @sprintf "Outputting x_l x_u = [%3.2f %3.2f]" x_l x_u)
        return (x_l, x_u)

    elseif f_l <= f_m <= f_u
        debug(LOGGER, "Case 2: Keep Moving Left")
        i = 1
        while f_l < f_m
            debug(LOGGER, @sprintf "x_u to x_m = [%3.2f, %3.2f]" x_l x_m)

            (x_u, x_m, f_u, f_m) = (x_m, x_l, f_m, f_l)
            x_l = x_l - (magnification^i) * abs(initial_step_length)
            f_l = f(x_l)

            i += 1
        end
        # f_l is now higher than f_m
        # This forms a valid bracketing interval
        info(LOGGER, @sprintf "Outputting x_l x_u = [%3.2f %3.2f]" x_l x_u)
        return (x_l, x_u)

    elseif f_l >= f_m <= f_u
        debug(LOGGER, "Case 3: Initial interval is bracket")
        info(LOGGER, @sprintf "Outputting x_l x_u = [%3.2f %3.2f]" x_l x_u)
        return (x_l, x_u)

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

    while !(F_c < F_a && F_c < F_b)
        debug(LOGGER, @sprintf "New loop iteration...")
        debug(LOGGER, @sprintf "F_a F_b F_c = %5.3f %5.3f %5.3f" F_a F_b F_c)
        a_current = a_next
        b_current = b_next
        c_current = c_next

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

        sleep(1)
        debug(LOGGER, @sprintf "a_next b_next c_next = %5.3f %5.3f %5.3f" a_next b_next c_next)
        debug(LOGGER, @sprintf "End of loop iteration...")
    end

    a_current = a_next
    b_current = b_next
    c_current = c_next

    debug(LOGGER, @sprintf "Exiting with a_current, b_current c_current = %5.3f %5.3f %5.3f" a_current b_current c_current)
    debug(LOGGER, @sprintf "Exiting with F_a F_b F_c = %5.3f %5.3f %5.3f" F_a F_b F_c)

    return (a_current, b_current)
end

end

