# Perform Fourier transform on MC, updating preallocated ηks matrix
function update_fourier!(mc::DimerMC)
    Lx, Ly = size(mc.spins)
    for (x, y) in Iterators.product(1:Lx, 1:Ly)
        disp = mc.spins[x,y] .- (x,y)
        mc.sks[x,y,:] .= [mod_equiv(disp, oriented_disps[i], mc) for i in 1:3]
    end
    fft!(mc.sks, (1, 2))
    mc.sks ./= length(mc.spins)

    for i in 1:3
        mc.ηks[:, :, i] .= getindex.(mc.ηs, i)
    end
    fft!(mc.ηks, (1, 2))
    mc.ηks ./= length(mc.ηs)
    return nothing
end

function validate_mc(mc::DimerMC)
    for I in eachindex(mc.spins)
        x,y = Tuple(I)
        @assert mod_equiv(mc.spins[mc.spins[x,y]...], (x,y), mc)
    end
end
