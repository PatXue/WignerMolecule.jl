# Perform Fourier transform on MC, updating preallocated spinks and ηks
# matrices, as well as calculating momentum-space correlations
function update_fourier!(mc::EtaMC)
    for i in 1:3
        mc.ηks[:, :, i] .= getindex.(mc.ηs, i)
    end
    fft!(mc.ηks, (1, 2))
    mc.ηks ./= length(mc.ηs)
    return nothing
end
