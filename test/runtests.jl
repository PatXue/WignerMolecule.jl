using WignerMolecule
using Random
using Test

paramtable = Dict(
    :re_Jpm => WignerParams(0, 0, 0, 0, 1.0, 0, 0, 0)
)
function init_dimermc(paramtype, ηi, ηj)
    wigparams = paramtable[paramtype]
    mc = DimerMC(; T=1.0, init_T=1.0, Q=0.5, wigparams, Lx=2, Ly=2)
    mc.ηs .= [ηi ηj; ηi ηj]
    mc.spins[1,2] = [2,1]
    mc.spins[2,1] = [1,2]
    addmonomer!((1,1), SpinVector(1,0,0), mc)
    addmonomer!((2,2), SpinVector(1,0,0), mc)
    return mc
end

@testset "WignerMolecule.jl" begin
    @testset "DimerMC J+- Re" begin
        mc = init_dimermc(:re_Jpm, SpinVector(1,0,0), SpinVector(1,0,0))
        @test validate_mc(mc)
    end
end
