# All energy units are in units of meV
struct WignerParams
    J_SS::Float64           # Spin-spin coupling
    J_EzEz_SS::Float64      # Spin-ηz coupling
    J_EzEz::Float64         # ηz coupling
    J_EAM_SS::ComplexF64    # Spin-η weird coupling (J+)
    J_EMEP_SS::ComplexF64   # Spin-η± coupling
    J_EMEM_SS::Float64      # Spin-η- coupling
    J_EMEP::ComplexF64      # η± coupling
    J_EMEM::Float64         # η- coupling
    H::Float64              # Bias field
end

# Default WignerParams values (for testing)
const default_params = WignerParams(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1e-3)

# Note: Using temperature in units of energy (k_B = 1)
struct WignerMC{AlgType} <: AbstractMC
    T::Float64     # Temperature
    params::WignerParams
    spins::PeriodicMatrix{SpinVector}
    ηs::PeriodicMatrix{SpinVector}

    spinks::Array{ComplexF64, 3}    # Fourier transformed spins
    ηks::Array{ComplexF64, 3}       # Fourier transformed ηs
    spink_corrs::Matrix{Float64}    # Correlations of Fourier spins
    ηk_corrs::Array{ComplexF64, 4}  # Correlations of Fourier ηs

    outdir::String # Output directory for local spin current plots
    savefreq::Int  # No. of sweeps between saving local spin current
end

function WignerMC{AlgType}(; T=1.0, wigparams=default_params, Lx=40, Ly=40,
    outdir=".", savefreq=0) where {AlgType}
    init_spins = fill(zeros(SpinVector), (Lx, Ly))
    init_ηs = fill(zeros(SpinVector), (Lx, Ly))
    return WignerMC{AlgType}(
        T, wigparams, init_spins, init_ηs,
        Array{ComplexF64}(undef, (Lx, Ly, 3)),
        Array{ComplexF64}(undef, (Lx, Ly, 3)),
        Matrix{Float64}(undef, (Lx, Ly)),
        Array{ComplexF64}(undef, (Lx, Ly, 3, 3)),
        outdir, savefreq
    )
end

function WignerMC{AlgType}(params::AbstractDict) where {AlgType}
    Lx, Ly = params[:Lx], params[:Ly]
    T = params[:T]
    wigparams = params[:wigparams]
    outdir = get(params, :outdir, ".")
    savefreq = get(params, :savefreq, 0)
    return WignerMC{AlgType}(; T, wigparams, Lx, Ly, outdir, savefreq)
end

function update_fourier!(mc::WignerMC)
    for i in 1:3
        mc.spinks[:, :, i] .= getindex.(mc.spins, i)
        mc.ηks[:, :, i] .= getindex.(mc.ηs, i)
    end
    fft!(mc.spinks, (1, 2))
    fft!(mc.ηks, (1, 2))
    mc.spink_corrs .= sum(abs2.(mc.spinks), dims=3)
    # ηk_iter = eachslice(mc.ηks, dims=(1, 2))
    # for I in eachindex(ηk_iter)
    #     η = ηk_iter[I]
    #     mc.ηk_corrs[I, :, :] .= η * η'
    # end
    return nothing
end

function Carlo.init!(mc::WignerMC, ctx::Carlo.MCContext, params::AbstractDict)
    init_type::Symbol = get(params, :init_type, :rand)
    if init_type == :const
        for I in eachindex(mc.spins)
            mc.spins[I] = SpinVector(0, 0, 1)
            mc.ηs[I] = SpinVector(0, 0, 1)
        end
    elseif init_type == :rand
        rand!(ctx.rng, mc.spins)
        rand!(ctx.rng, mc.ηs)
    elseif init_type == :afm_fe
        init_afm_fe!(mc.spins, mc.ηs)
    end
    update_fourier!(mc)
    return nothing
end

function Carlo.sweep!(mc::WignerMC, ctx::Carlo.MCContext)
    Carlo.sweep!(mc, ctx.rng)
end

const ω::ComplexF64 = exp(im * 2π/3)
# Calculate the energy at a lattice site (x, y) if it had spin s and
# pseudospin η
function energy(mc::WignerMC, s::SpinVector, η::SpinVector, x, y)
    # Nearest neighbor lattice positions
    nns = ((x+1, y), (x+1, y-1), (x, y-1), (x-1, y), (x-1, y+1), (x, y+1))
    # Coupling energies
    J_SS = mc.params.J_SS
    J_EzEz_SS = mc.params.J_EzEz_SS
    J_EzEz = mc.params.J_EzEz
    J_EAM_SS = mc.params.J_EAM_SS
    J_EMEP_SS = mc.params.J_EMEP_SS
    J_EMEM_SS = mc.params.J_EMEM_SS
    J_EMEP = mc.params.J_EMEP
    J_EMEM = mc.params.J_EMEM

    E = 0.0
    for j in eachindex(nns)
        ν = ω^(j-1)
        sj = mc.spins[nns[j]...]
        ηj = mc.ηs[nns[j]...]

        # η raising and lowering operators
        η_m = η[1] + 1.0im*η[2]
        ηj_p = ηj[1] - 1.0im*ηj[2]
        ηj_m = ηj[1] + 1.0im*ηj[2]

        E_spin = 0.0 + 0.0im
        E_η = 0.0 + 0.0im

        # η-only energy
        E_η +=   J_EzEz *     η[3] * ηj[3]
        E_η += 2*J_EMEP *     η_m * ηj_p
        E_η += 2*J_EMEM * ν * η_m * ηj_m

        # η-S energy
        E_spin +=   J_EzEz_SS *     η[3] * ηj[3]
        E_spin += 2*J_EMEP_SS *     η_m * ηj_p
        E_spin += 2*J_EMEM_SS * ν * η_m * ηj_m
        E_spin += J_SS
        E_spin += 2*J_EAM_SS * (η_m/ν + ηj_p*ν)

        E += real(E_spin * (s⋅sj) + E_η)
    end
    E -= mc.params.H * s[3]

    return E
end

# Calculate the energy contribution of a site (x, y), considering only half of
# its bonds (avoids double counting when calculating total energy)
function half_energy(mc::WignerMC, x, y)
    # Nearest neighbor lattice positions
    nns = ((x+1, y), (x+1, y-1), (x, y-1))
    # Coupling energies
    J_SS = mc.params.J_SS
    J_EzEz_SS = mc.params.J_EzEz_SS
    J_EzEz = mc.params.J_EzEz
    J_EAM_SS = mc.params.J_EAM_SS
    J_EMEP_SS = mc.params.J_EMEP_SS
    J_EMEM_SS = mc.params.J_EMEM_SS
    J_EMEP = mc.params.J_EMEP
    J_EMEM = mc.params.J_EMEM

    s = mc.spins[x, y]
    η = mc.ηs[x, y]
    E = 0.0
    for j in eachindex(nns)
        ν = ω^(j-1)
        sj = mc.spins[nns[j]...]
        ηj = mc.ηs[nns[j]...]

        # η raising and lowering operators
        η_m = η[1] + 1.0im*η[2]
        ηj_p = ηj[1] - 1.0im*ηj[2]
        ηj_m = ηj[1] + 1.0im*ηj[2]

        E_spin = 0.0 + 0.0im
        E_η = 0.0 + 0.0im

        # η-only energy
        E_η +=   J_EzEz *     η[3] * ηj[3]
        E_η += 2*J_EMEP *     η_m * ηj_p
        E_η += 2*J_EMEM * ν * η_m * ηj_m

        # η-S energy
        E_spin +=   J_EzEz_SS *     η[3] * ηj[3]
        E_spin += 2*J_EMEP_SS *     η_m * ηj_p
        E_spin += 2*J_EMEM_SS * ν * η_m * ηj_m
        E_spin += J_SS
        E_spin += 2*J_EAM_SS * (η_m/ν + ηj_p*ν)

        E += real(E_spin * (s⋅sj) + E_η)
    end
    E -= mc.params.H * s[3]

    return E
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

function Carlo.measure!(mc::WignerMC, ctx::Carlo.MCContext)
    Lx, Ly = size(mc.spins)
    N = Lx * Ly
    # Magnetization per lattice site
    mag = norm(sum(mc.spins)) / N
    measure!(ctx, :Mag, mag)
    measure!(ctx, :Mag2, mag^2)
    measure!(ctx, :Mag4, mag^4)

    η_sum = sum(mc.ηs)
    measure!(ctx, :ηz, η_sum[3] / N)
    measure!(ctx, :ηxy, sqrt(η_sum[1]^2 + η_sum[2]^2) / N)

    # Energy per lattice site
    energy = 0.0
    for y in 1:Ly
        for x in 1:Lx
            energy += half_energy(mc, x, y)
        end
    end
    energy /= N
    measure!(ctx, :Energy, energy)
    measure!(ctx, :Energy2, energy^2)

    update_fourier!(mc)
    measure!(ctx, :spink_corrs, mc.spink_corrs)
    # measure!(ctx, :etak_corrs, mc.ηk_corrs)

    if is_save_sweep(mc, ctx)
        save_spin(mc, ctx)
    end

    return nothing
end

function Carlo.register_evaluables(::Type{WignerMC{AlgType}}, eval::AbstractEvaluator,
                                   params::AbstractDict) where {AlgType}
    T = params[:T]
    N = params[:Lx] * params[:Ly]
    evaluate!(eval, :χ, (:Mag, :Mag2)) do mag, mag2
        return N / T * (mag2 - mag^2)
    end

    evaluate!(eval, :HeatCap, (:Energy2, :Energy)) do E2, E
        return N * (E2 - E^2) / T^2
    end

    return nothing
end

function Carlo.write_checkpoint(mc::WignerMC, out::HDF5.Group)
    out["spins"] = mc.spins
    out["etas"] = mc.ηs
    return nothing
end
function Carlo.read_checkpoint!(mc::WignerMC, in::HDF5.Group)
    raw_spins = read(in, "spins")
    raw_ηs = read(in, "etas")
    mc.spins .= map(v -> SVector(v[:data][1], v[:data][2], v[:data][3]), raw_spins)
    mc.ηs .= map(v -> SVector(v[:data][1], v[:data][2], v[:data][3]), raw_ηs)
    return nothing
end
