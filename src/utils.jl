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

# Calculate temperature during thermalization
function calc_temp(mc, ctx::Carlo.MCContext)
    if is_thermalized(ctx)
        return mc.T
    else
        return mc.init_T + (mc.T - mc.init_T) * ctx.sweeps/ctx.thermalization_sweeps
    end
end
# Calculate bias field during thermalization
function calc_B(mc::WignerMC, ctx::Carlo.MCContext)
    if is_thermalized(ctx)
        return mc.B
    else
        ΔB = (mc.B - mc.init_B) * ctx.sweeps/ctx.thermalization_sweeps * 2
        return mc.init_B + sign(ΔB) * min(abs(ΔB), abs(mc.B - mc.init_B))
    end
end

# Indexing functions for Γ and M point correlators
Γ(corrs) = corrs[1,1,:]
function M(corrs)
    Lx, _ = size(corrs)
    return corrs[div(Lx, 2)+1, 1, :]
end
function M2(corrs)
    _, Ly = size(corrs)
    return corrs[1, div(Ly,2)+1, :]
end
function M3(corrs)
    Lx, Ly = size(corrs)
    return corrs[div(Lx,2)+1, div(Ly,2)+1, :]
end
function half_M(corrs)
    Lx, _ = size(corrs)
    return corrs[div(Lx, 4)+1, 1, :]
end
function half_M2(corrs)
    _, Ly = size(corrs)
    return corrs[1, div(Ly, 4)+1, :]
end
function half_M3(corrs)
    Lx, Ly = size(corrs)
    return corrs[div(Lx, 4)+1, div(Ly, 4)+1, :]
end
function K(corrs)
    Lx, _ = size(corrs)
    n = div(Lx, 3)
    return corrs[2n+1,n+1,:]
end
function part_K(corrs) # 3/4 of the K point, not sure why
    Lx, _ = size(corrs)
    n = div(Lx, 4)
    return corrs[2n+1,n+1,:]
end
function part_K2(corrs)
    Lx, _ = size(corrs)
    n = div(Lx, 4)
    return corrs[3n+1,n+1,:]
end
function part_K3(corrs)
    Lx, _ = size(corrs)
    n = div(Lx, 4)
    return corrs[3n+1,2n+1,:]
end
const corr_posns = (Γ, M, M2, M3, half_M, half_M2, half_M3, K, part_K, part_K2, part_K3)

# Calculate sum of chirality for each triangular plaquette
function chirality(spins)
    Q = 0.0
    for I in eachindex(spins)
        x, y = Tuple(I)
        s = spins[I]
        Q += s ⋅ (spins[x+1,y] × spins[x,y+1])
    end
    return Q
end
# Calculate chirality for a single site
function chirality(spins, x, y)
    return spins[x,y] ⋅ (spins[x+1,y] × spins[x,y+1])
end

norm2(v) = sum(abs2.(v))

# Randomize spins in MC
function randomize!(mc::WignerMC)
    rand!(mc.spins)
    rand!(mc.ηs)
    update_fourier!(mc)
end
