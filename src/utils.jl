const SpinVector = SVector{3, Float64}
function Random.rand(rng::AbstractRNG, ::Random.SamplerType{SpinVector})
    ϕ = 2π * rand(rng)
    θ = acos(2 * rand(rng) - 1)
    return SpinVector(cos(ϕ)sin(θ), sin(ϕ)sin(θ), cos(θ))
end

# Sum spins of position (x, y)'s nearest neighbors
nn_sum(A::AbstractArray, x, y) = A[x+1, y] + A[x, y+1] + A[x-1, y] + A[x, y-1]
# Sum spins of (x, y)'s next nearest NE-SW neighbors
nnna_sum(A::AbstractArray, x, y) = A[x+1, y+1] + A[x-1, y-1]
# Sum spins of (x, y)'s next nearest NW-SE neighbors
nnnb_sum(A::AbstractArray, x, y) = A[x+1, y-1] + A[x-1, y+1]
