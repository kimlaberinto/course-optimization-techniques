function NDRosenbrock(N::Integer, x::Array{T}) where T <: Real
    @assert length(x) == N
    @assert N >= 2
    result = 0
    for i in 1:(N-1)
        result += 100*(x[i+1] - x[i]^2)^2 + (1 - x[i])^2
    end
    return result
end