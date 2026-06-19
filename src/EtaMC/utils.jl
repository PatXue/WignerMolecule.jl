# Perform Fourier transform on MC, updating preallocated spinks and ηks
# matrices, as well as calculating momentum-space correlations
function update_fourier!(mc::EtaMC)
    for i in 1:3
        mc.spinks[:, :, i] .= getindex.(mc.spins, i)
    end
    fft!(mc.spinks, (1, 2))
    mc.spinks ./= length(mc.spins)

    for I in eachindex(IndexCartesian(), mc.spins)
        x, y = Tuple(I)
        mc.chis[I] = chirality(mc.spins, x, y)
    end
    fft!(mc.chis)
    mc.chis ./= length(mc.spins)
    return nothing
end
