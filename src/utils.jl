function rand_spin(rng::AbstractRNG = default_rng())
    ϕ = 2π * rand(rng)
    θ = acos(2 * rand(rng) - 1)
    return (cos(ϕ)sin(θ), sin(ϕ)sin(θ), cos(θ))
end

function init_eag!(spins::AbstractMatrix{SpinVector})
    for I in eachindex(IndexCartesian(), spins)
        x, y = Tuple(I)
        spin_sign = mod(x+y, 4) < 2 ? 1.0 : -1.0
        spins[I] = spin_sign * SVector(1.0, 0.0, 0.0)
    end
end

function init_orth!(spins::AbstractMatrix{SpinVector})
    for I in eachindex(IndexCartesian(), spins)
        x, y = Tuple(I)
        θ = π/2 * (x + y)
        spins[I] = SVector(cos(θ), sin(θ), 0.0)
    end
end

function init_afm_fe!(spins::AbstractMatrix{SpinVector}, ηs::AbstractMatrix{SpinVector})
    for I in eachindex(IndexCartesian(), spins)
        x, _ = Tuple(I)
        spins[I] = SVector(0, 0, (-1)^x)
        ηs[I] = SVector(cos(π/3), sin(π/3), 0)
    end
end
