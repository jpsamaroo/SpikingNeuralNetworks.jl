@with_kw struct HHParameter{VF}
    Cm::VF = 1uF * cm^(-2) * 20000um^2
    gl::VF = 5e-5siemens * cm^(-2) * 20000um^2
    El::VF = -65mV
    Ek::VF = -90mV
    En::VF = 50mV
    gn::VF = 100msiemens * cm^(-2) * 20000um^2
    gk::VF = 30msiemens * cm^(-2) * 20000um^2
    Vt::VF = -63mV
    τe::VF = 5ms
    τi::VF = 10ms
    Ee::VF = 0mV
    Ei::VF = -80mV
end

@with_kw mutable struct HH{VF,VB}
    param::HHParameter = HHParameter()
    N::SNNInt = 100
    v::VF = param.El .+ 5(randn(eltype(VF), N) .- 1)
    m::VF = zeros(eltype(VF), N)
    n::VF = zeros(eltype(VF), N)
    h::VF = ones(eltype(VF), N)
    ge::VF = (1.5randn(eltype(VF), N) .+ 4) .* 10nS
    gi::VF = (12randn(eltype(VF), N) .+ 20) .* 10nS
    fire::VB = zeros(eltype(VB), N)
    I::VF = zeros(eltype(VF), N)
    records::Dict = Dict()
end
HH(x;kwargs...) = HH{Vector{SNNFloat},Vector{Bool}}(;kwargs...)

function integrate!(p::HH, param::HHParameter, dt::SNNFloat)
    @unpack N, v, m, n, h, ge, gi, fire, I = p
    @unpack Cm, gl, El, Ek, En, gn, gk, Vt, τe, τi, Ee, Ei = param
    @inbounds for i = 1:N
        m[i] += dt * (0.32f0 * (13f0 - v[i] + Vt) / (exp((13f0 - v[i] + Vt) / 4f0) - 1f0) * (1f0 - m[i]) -
        0.28f0 * (v[i] - Vt - 40f0) / (exp((v[i] - Vt - 40f0) / 5f0) - 1f0) * m[i])
        n[i] += dt * (0.032f0 * (15f0 - v[i] + Vt) / (exp((15f0 - v[i] + Vt) / 5f0) - 1f0) * (1f0 - n[i]) -
        0.5f0 * exp((10f0 - v[i] + Vt) / 40f0) * n[i])
        h[i] += dt * (0.128f0 * exp((17f0 - v[i] + Vt) / 18f0) * (1f0 - h[i]) -
        4f0 / (1f0 + exp((40f0 - v[i] + Vt) / 5f0)) * h[i])
        v[i] += dt / Cm * ( I[i] + gl * (El - v[i]) + ge[i] * (Ee - v[i]) + gi[i] * (Ei - v[i]) +
        gn * m[i]^3 * h[i] * (En - v[i]) + gk * n[i]^4 * (Ek - v[i]) )
        ge[i] += dt * -ge[i] / τe
        gi[i] += dt * -gi[i] / τi
    end
    @inbounds for i = 1:N
        fire[i] = v[i] > -20f0
    end
end
