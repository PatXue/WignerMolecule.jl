# Initialize the DimerMC spins/etas
const unitcellspins = [(0,1) (0,-1); (1,-1) (-1,1)]
function init_vbs_s!(mc::DimerMC)
    for I in eachindex(IndexCartesian(), mc.spins)
        x, y = Tuple(I)
        # Shift x,y to be inside unit cell
        while y > 2
            x += 1
            y -= 2
        end
        x = mod1(x, 2)
        mc.spins[I] = SVector((Tuple(I) .+ unitcellspins[x,y])...)
    end
end

const unitcelletas = [(-0.81,0.31,0) (0.52,-0.61,0); (0.49,0.61,0) (-0.78,-0.34,0)]
function init_vbs_eta!(mc::DimerMC)
    for I in eachindex(IndexCartesian(), mc.ηs)
        x, y = Tuple(I)
        # Shift x,y to be inside unit cell
        while y > 2
            x += 1
            y -= 2
        end
        x = mod1(x, 2)
        eta = unitcelletas[x,y] ./ norm(unitcelletas[x,y])
        mc.ηs[I] = SVector(eta...)
    end
end

function Carlo.init!(mc::DimerMC, ctx::Carlo.MCContext, params::AbstractDict)
    init_type = get(params, :init_type, :vbs)
    empty!(mc.monomers)
    if init_type == :vbs
        init_vbs_s!(mc)
        init_vbs_eta!(mc)
    elseif init_type == :rand
        rand!(ctx.rng, mc.monospins)
        rand!(ctx.rng, mc.ηs)
        for I in eachindex(mc.spins)
            addmonomer!(convert(SVector, I), mc.monospins[I], mc)
        end
    elseif init_type == :fm
        for I in eachindex(mc.spins)
            addmonomer!(convert(SVector, I), SVector(0,0,1), mc)
            mc.ηs[I] = SVector(0,0,1)
        end
    elseif init_type == :stripe
        init_stripe!(mc.monospins, mc.ηs)
        for I in eachindex(mc.spins)
            addmonomer!(convert(SVector, I), mc.monospins[I], mc)
        end
    end
    update_fourier!(mc)
    return nothing
end
