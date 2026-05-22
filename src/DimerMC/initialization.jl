# Initialize the DimerMC spins/etas
const unitcell = [a2 a6; a5 a3]
function Carlo.init!(mc::DimerMC, ctx::Carlo.MCContext, params::AbstractDict)
    for I in eachindex(IndexCartesian(), mc.spins)
        x, y = Tuple(I)
        # Shift x,y to be inside unit cell
        while y > 2
            x += 1
            y -= 2
        end
        x = (x % 2) + 1

        mc.spins[I] = unitcell[x,y]
    end

    init_type::Symbol = get(params, :init_type, :rand)
    rand!(ctx.rng, mc.ηs)
    if init_type == :const
        for I in eachindex(mc.ηs)
            mc.ηs[I] = SpinVector(0, 0, 1)
        end
    elseif init_type == :afm_fe_eta
        init_afm_fe_eta!(mc.ηs)
    elseif init_type == :stripe_eta
        init_stripe_eta!(mc.ηs)
    end

    update_fourier!(mc)
    return nothing
end
