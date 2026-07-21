using WignerMolecule
using Random
using Test

@testset "WignerMolecule.jl" begin
    @testset "DimerMC J+- Re" begin
        wigparams = WignerParams(0, 0, 0, 0, 1.0, 0, 0, 0)
        mc = DimerMC(; T=1.0, init_T=1.0, Q=0.5, wigparams, Lx=2, Ly=2)
        mc.ηs .= fill(SpinVector(1,0,0), 2, 2)
        mc.spins[1,2] = [2,1]
        mc.spins[2,1] = [1,2]
        addmonomer!((1,1), SpinVector(1,0,0), mc)
        addmonomer!((2,2), SpinVector(1,0,0), mc)
        @test validate_mc(mc)
    end
end
