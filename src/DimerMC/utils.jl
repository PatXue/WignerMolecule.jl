# Perform Fourier transform on MC, updating preallocated ηks matrix
function update_fourier!(mc::DimerMC)
    for i in 1:3
        mc.ηks[:, :, i] .= getindex.(mc.ηs, i)
    end
    fft!(mc.ηks, (1, 2))
    mc.ηks ./= length(mc.ηs)
    return nothing
end
