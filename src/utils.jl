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

# Perform Fourier transform on MC, updating preallocated spinks and ηks
# matrices, as well as calculating momentum-space correlations
function update_fourier!(mc::WignerMC)
    for i in 1:3
        mc.spinks[:, :, i] .= getindex.(mc.spins, i)
        mc.ηks[:, :, i] .= getindex.(mc.ηs, i)
    end
    fft!(mc.spinks, (1, 2))
    mc.spinks ./= length(mc.spins)
    fft!(mc.ηks, (1, 2))
    mc.ηks ./= length(mc.ηs)
    mc.spink_corrs .= sum(abs2.(mc.spinks), dims=3)
    ηk_iter = eachslice(mc.ηks, dims=(1, 2))
    for I in eachindex(ηk_iter)
        η = ηk_iter[I]
        mc.ηk_corrs[I, :, :] .= η .* η'
    end
    return nothing
end

# Returns if spins should be saved on this sweep (assuming thermalized)
function is_save_sweep(mc::WignerMC, ctx::Carlo.MCContext)
    measure_sweeps = ctx.sweeps - ctx.thermalization_sweeps
    return mc.savefreq > 0 && measure_sweeps % mc.savefreq == 0
end

# Save spins to JLD2 file
function save_spin(mc::WignerMC, ctx::Carlo.MCContext)
    jldopen("$(mc.outdir)/spins.jld2", "a+") do file
        sweep_name = "sweep$(ctx.sweeps - ctx.thermalization_sweeps)"
        file[sweep_name * "-spins"] = Matrix(mc.spins)
        file[sweep_name * "-etas"] = Matrix(mc.ηs)
    end
end

function calc_temp(mc::WignerMC, ctx::Carlo.MCContext)
    if is_thermalized(ctx)
        return mc.T
    else
        return mc.init_T + (mc.T - mc.init_T) * ctx.sweeps/ctx.thermalization_sweeps
    end
end

function calc_B(mc::WignerMC, ctx::Carlo.MCContext)
    if is_thermalized(ctx)
        return mc.B
    else
        return mc.init_B + (mc.B - mc.init_B) * ctx.sweeps/ctx.thermalization_sweeps
    end
end
