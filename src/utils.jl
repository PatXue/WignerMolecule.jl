"""
    rand_spin(rng::AbstractRNG = default_rng())

Generate a uniformly random spin vector (normalized to 1)
"""
function rand_spin(rng::AbstractRNG = default_rng())
    ϕ = 2π * rand(rng)
    θ = acos(2 * rand(rng) - 1)
    return (cos(ϕ)sin(θ), sin(ϕ)sin(θ), cos(θ))
end

"""
    spins_iter(mc::WignerMC)

Returns iterable over spin-η pairs in `mc`
"""
function spins_iter(mc::WignerMC)
    zip(eachslice(mc.spins, dims=(2, 3)), eachslice(mc.ηs, dims=(2, 3)))
end

sites_iter(mc::WignerMC) = Iterators.product(axes(mc.spins, 2), axes(mc.spins, 3))

function init_orth!(mc::WignerMC)
    for (x, y) in sites_iter(mc)
        θ = π/2 * (x + y)
        mc.spins[:, x, y] .= (cos(θ), sin(θ), 0.0)
    end
end

function init_afm_fe!(mc::WignerMC)
    for (x, y) in sites_iter(mc)
        mc.spins[:, x, y] .= (0, 0, (-1)^x)
        mc.ηs[:, x, y] .= (cos(pi/2), -sin(pi/2), 0)
    end
end
