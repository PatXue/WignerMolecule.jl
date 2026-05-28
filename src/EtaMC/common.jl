# Initialize the WignerMC spins
function Carlo.init!(mc::EtaMC, ctx::Carlo.MCContext, params::AbstractDict)
    init_type::Symbol = get(params, :init_type, :rand)

    rand!(ctx.rng, mc.spins)
    if init_type == :fm
        for I in eachindex(mc.spins)
            mc.spins[I] = SpinVector(0, 0, 1)
        end
    elseif init_type == :fe
        for I in eachindex(mc.spins)
            mc.spins[I] = SpinVector(0, 1, 0)
        end
    elseif init_type == :stripe
        for I in eachindex(mc.spins)
            _, y = Tuple(I)
            mc.spins[I] = SpinVector((-1)^y, 0, 0)
        end
    end

