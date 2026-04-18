# Sampling functions for high-T expectations
module HighTemp

# Randomize spins in MC, updates Fourier transform if flag passed
function randomize!(mc::WignerMC; up_fourier=false)
    rand!(mc.spins)
    rand!(mc.ηs)
    if up_fourier
        update_fourier!(mc)
    end
end

end
