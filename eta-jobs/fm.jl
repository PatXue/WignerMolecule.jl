import Pkg
Pkg.activate("..")

using WignerMolecule

tm = TaskMaker()
jobname = "fm"
tm.wigparams = EtaParams(1.5, 0.5)

Ls = [20]
Ts = 0.1:0.1:1.0
for (T, L) in Iterators.product(Ts, Ls)
    tm.T = T
    tm.Lx = tm.Ly = L
    task(tm)
end
