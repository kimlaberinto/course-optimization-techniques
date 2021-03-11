module objectivefunctionModule

using ForwardDiff

export NDRosenbrock

export analyticGradientNDRosenbrock
export analyticHessianNDRosenbrock

export autodiffGradientNDRosenbrock
export autodiffHessianNDRosenbrock

function NDRosenbrock(N::Integer, x::Array{T}) where T <: Real
    @assert length(x) == N
    @assert N >= 2
    result = 0
    for i in 1:(N-1)
        result += 100*(x[i+1] - x[i]^2)^2 + (1 - x[i])^2
    end
    return result
end

typeATerm(x_i) = 2*(1-x_i)*(-1)
typeBTerm(x_i, x_next) = 100*2*(x_next - x_i^2)*(-2*x_i)
typeCTerm(x_prev, x_i) = 100*2*(x_i - x_prev^2)
function analyticGradientNDRosenbrock(N::Integer, x::Array{T}) where T <: Real
    gradresult = zeros(eltype(x), N)

    for i in 1:N
        if i == 1 #First
            gradresult[i] = 
                typeATerm(x[i]) + typeBTerm(x[i], x[i+1])
        elseif i == N #Last
            gradresult[i] =
                typeCTerm(x[i-1], x[i])
        else #Intermediary terms
            gradresult[i] = 
                typeATerm(x[i]) + typeBTerm(x[i], x[i+1]) + typeCTerm(x[i-1], x[i])
        end
    end

    return gradresult
end

function analyticHessianNDRosenbrock(N::Integer, x::Array{T}) where T <: Real
    error("Currently undefined. Please use autodiffHessianNDRosenbrock.")
end

function autodiffGradientNDRosenbrock(N::Integer, x::Array{T}) where T <: Real
    return ForwardDiff.gradient(in -> NDRosenbrock(N, in), x)
end

function autodiffHessianNDRosenbrock(N::Integer, x::Array{T}) where T <: Real
    return ForwardDiff.hessian(in -> NDRosenbrock(N, in), x)
end

end

