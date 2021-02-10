module A1Module

using Memento #Invenia Module for logging
using Printf

export SwannsBracketingMethod

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

end

