# Initialize the DimerMC spins/etas
function Carlo.init!(mc::DimerMC, ctx::Carlo.MCContext, params::AbstractDict)
    init_type::Symbol = get(params, :init_type, :rand)
    rand!(ctx.rng, mc.ηs)
    if init_type == :const
        for I in eachindex(mc.spins)
            mc.spins[I] = SpinVector(0, 0, 1)
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
