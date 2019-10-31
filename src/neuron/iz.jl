@with_kw struct IZParameter
    a::SNNFloat = 0.01
    b::SNNFloat = 0.2
    c::SNNFloat = -65
    d::SNNFloat = 2
end

@with_kw mutable struct IZ{VF,VB}
    param::IZParameter = IZParameter()
    N::SNNInt = 100
    v::VF = fill(-65.0, N)
    u::VF = param.b * v
    fire::VB = zeros(eltype(VB), N)
    I::VF = zeros(eltype(VF), N)
    records::Dict = Dict()
end
IZ(x;kwargs...) = IZ{Vector{SNNFloat},Vector{Bool}}(;kwargs...)

function integrate!(p::IZ, param::IZParameter, dt::SNNFloat)
    @unpack N, v, u, fire, I = p
    @unpack a, b, c, d = param
    @inbounds for i = 1:N
        v[i] += 0.5f0dt * (0.04f0v[i]^2 + 5f0v[i] + 140f0 - u[i] + I[i])
        v[i] += 0.5f0dt * (0.04f0v[i]^2 + 5f0v[i] + 140f0 - u[i] + I[i])
        u[i] += dt * (a * (b * v[i] - u[i]))
    end
    @inbounds for i = 1:N
        fire[i] = v[i] > 30f0
        v[i] = ifelse(fire[i], c, v[i])
        u[i] += ifelse(fire[i], d, 0f0)
    end
end
