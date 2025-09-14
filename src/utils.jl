const SpinVector = SVector{3, Float64}
function Random.rand(rng::AbstractRNG, ::Random.SamplerType{SpinVector})
    ϕ = 2π * rand(rng)
    θ = acos(2 * rand(rng) - 1)
    return SpinVector(cos(ϕ)sin(θ), sin(ϕ)sin(θ), cos(θ))
end
