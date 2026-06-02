# Perform Fourier transform on MC, updating preallocated ηks matrix
function update_fourier!(mc::DimerMC)
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
