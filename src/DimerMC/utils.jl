# Perform Fourier transform on MC, updating preallocated ηks matrix
function update_fourier!(mc::DimerMC)
    Lx, Ly = size(mc.spins)
    for (x, y) in Iterators.product(1:Lx, 1:Ly)
        disp = mc.spins[x,y] .- (x,y)
        # Map bonds along a1,a2,a3 to 1,2,3 and -1,-2,-3 for -a1,-a2,-a3
        for i in 1:3
            if mod_equiv(disp, oriented_disps[i], mc)
                mc.sks[x,y] = i
            elseif mod_equiv(disp, -oriented_disps[i], mc)
                mc.sks[x,y] = -i
            end
        end
    end
    fft!(mc.sks)
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
