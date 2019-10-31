abstract type AbstractIFParameter end
@with_kw struct IFParameter{VF} <: AbstractIFParameter
    τm::VF = 20ms
    τe::VF = 5ms
    τi::VF = 10ms
    Vt::VF = -50mV
    Vr::VF = -60mV
    El::VF = Vr
end

abstract type AbstractIF end
@with_kw mutable struct IF{VF,VB} <: AbstractIF
    param::IFParameter = IFParameter()
    N::SNNInt = 100
    v::VF = param.Vr .+ rand(eltype(VF), N) .* (param.Vt - param.Vr)
    ge::VF = zeros(eltype(VF), N)
    gi::VF = zeros(eltype(VF), N)
    fire::VB = zeros(eltype(VB), N)
    I::VF = zeros(eltype(VF), N)
    records::Dict = Dict()
end
IF(x;kwargs...) = IF{Vector{SNNFloat},Vector{Bool}}(;kwargs...)

function integrate!(p::IF, param::IFParameter, dt::SNNFloat)
    @unpack N, v, ge, gi, fire, I = p
    @unpack τm, τe, τi, Vt, Vr, El = param
    @inbounds for i = 1:N
        v[i] += dt * (ge[i] + gi[i] - (v[i] - El) + I[i]) / τm
        ge[i] += dt * -ge[i] / τe
        gi[i] += dt * -gi[i] / τi
    end
    @inbounds for i = 1:N
        fire[i] = v[i] > Vt
        v[i] = ifelse(fire[i], Vr, v[i])
    end
end
