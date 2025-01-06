module RandomGeneratorsModule

export exponential, pareto, erlang

using Random

function exponential(λ::Float64)
    return rand(Exponential(1.0 / λ))
end

function pareto(α::Float64, scale::Float64=1.0)
    return (rand(Pareto(α)) + 1) * scale
end

function erlang(k::Int, λ::Float64)
    return sum(rand(Exponential(1.0 / λ)) for _ in 1:k)
end

end # module
