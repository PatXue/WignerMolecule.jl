using WignerMolecule
const wm = WignerMolecule
using Random
using Test

paramtable = Dict(
    :re_Jpmss => WignerParams(0, 0, 0, 0, 1.0, 0, 0, 0)
)
function init_dimermc(paramtype, η)
    wigparams = paramtable[paramtype]
    mc = DimerMC(; T=1.0, init_T=1.0, Q=0.5, wigparams, Lx=2, Ly=2)
    fill!(mc.ηs, η)
    mc.spins[1,2] = [2,1]
    mc.spins[2,1] = [1,2]
    addmonomer!((1,1), SpinVector(1,0,0), mc)
    addmonomer!((2,2), SpinVector(1,0,0), mc)
    return mc
end

@testset "WignerMolecule.jl" begin
    @testset "DimerMC Jpmss Re" begin
        mc = init_dimermc(:re_Jpmss, SpinVector(1,0,0))
        @test validate_mc(mc)
        d = wm.Dimer([1,2], [2,1])
        s = SpinVector(1,0,0)
        dimerE = -3/4 * 2 * 1/4
        siteE = 2 * 2 * 1/4 * 1/4

        @test wm.dimer_energy_s(mc, d) == dimerE
        @test wm.ssfactor(mc, wm.Dimer([1,1], [2,2])) == 2 * 1/4
        @test wm.bond_energy_s(mc, wm.Dimer([1,1], [2,2]), 1/4) == siteE / 2
        @test wm.bond_energy(mc, wm.Dimer([1,1], [2,2])) == siteE / 2
        @test wm.site_energy_s(mc, [1,1], s) == siteE
        @test wm.shift_energy_s(mc, d, [1,1], s) == dimerE + siteE
        @test wm.pair_energy_s(mc, [1,1], [2,2], s, s) == siteE / 2
        @test wm.half_energy(mc, [1,1]) == siteE / 2
        @test wm.half_energy(mc, [1,2]) == dimerE
        @test wm.half_energy(mc, [2,1]) == dimerE
        @test wm.half_energy(mc, [2,2]) == siteE / 2
        @test total_energy(mc) == 2*dimerE + siteE
    end
end
