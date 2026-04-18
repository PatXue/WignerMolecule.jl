# Uniform random 3D unit vector
function Random.rand(rng::AbstractRNG, ::Random.SamplerType{SpinVector})
    ϕ = 2π * rand(rng)
    θ = acos(2 * rand(rng) - 1)
    return SpinVector(cos(ϕ)sin(θ), sin(ϕ)sin(θ), cos(θ))
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

function init_afm_fe_s!(spins::AbstractMatrix{SpinVector})
    for I in eachindex(IndexCartesian(), spins)
        x, _ = Tuple(I)
        spins[I] = SVector(0, 0, (-1)^x)
    end
end

function init_afm_fe_eta!(ηs::AbstractMatrix{SpinVector})
    for I in eachindex(IndexCartesian(), ηs)
        ηs[I] = SVector(cos(π/3), -sin(π/3), 0)
    end
end

function init_afm_fe!(spins::AbstractMatrix{SpinVector}, ηs::AbstractMatrix{SpinVector})
    init_afm_fe_s!(spins)
    init_afm_fe_eta!(ηs)
end

function init_stripe_s!(spins::AbstractMatrix{SpinVector})
    for I in eachindex(IndexCartesian(), spins)
        x, _ = Tuple(I)
        spins[I] = SVector(0, 0, (-1)^(x ÷ 2))
    end
end

function init_stripe_eta!(ηs::AbstractMatrix{SpinVector})
    for I in eachindex(IndexCartesian(), ηs)
        x, _ = Tuple(I)
        ηs[I] = (-1)^x .* SVector(cos(π/3), -sin(π/3), 0)
    end
end

function init_stripe!(spins::AbstractMatrix{SpinVector}, ηs::AbstractMatrix{SpinVector})
    init_stripe_s!(spins)
    init_stripe_eta!(ηs)
end

function init_afm_afe_s!(spins::AbstractMatrix{SpinVector})
    for I in eachindex(IndexCartesian(), spins)
        x, y = Tuple(I)
        spins[I] = SVector(0, 0, (-1)^(x + div(y,2)))
    end
end

function init_afm_afe_eta!(ηs::AbstractMatrix{SpinVector})
    for I in eachindex(IndexCartesian(), ηs)
        _, y = Tuple(I)
        ηs[I] = SVector(0, (-1)^(y+1), 0)
    end
end

function init_afm_afe!(spins::AbstractMatrix{SpinVector}, ηs::AbstractMatrix{SpinVector})
    init_afm_afe_s!(spins)
    init_afm_afe_eta!(ηs)
end

function init_skyrm!(mc::WignerMC, xi; etatype=:const)
    if etatype == :const
        fill!(mc.ηs, SpinVector(0, 0, 1))
    elseif etatype == :stripe
        init_stripe_eta!(mc.ηs)
    else
        rand!(mc.ηs)
    end

    Lx, Ly = size(mc.spins)
    X, Y = Lx/2, Ly/2
    for (x, y) in Iterators.product(1:Lx, 1:Ly)
        r = sqrt((x-X)^2 + (y-Y)^2)
        if r > xi
            mc.spins[x, y] = SpinVector(0, 0, 1)
        elseif r == 0
            mc.spins[x, y] = SpinVector(0, 0, -1)
        else
            mc.spins[x, y] = SpinVector(sin(π*r/xi) * ((x-X)/r), sin(π*r/xi) * ((y-Y)/r), -cos(π*r/xi))
        end
    end
end
