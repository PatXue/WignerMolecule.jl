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
    if init_type == :vbs
        init_vbs_s!(mc)
        init_vbs_eta!(mc)
        empty!(mc.monomers)
    elseif init_type == :rand
        for I in eachindex(mc.spins)
            mc.spins[I] = I
        end
        rand!(ctx.rng, mc.monospins)
        rand!(ctx.rng, mc.ηs)
    end
    update_fourier!(mc)
    return nothing
end
