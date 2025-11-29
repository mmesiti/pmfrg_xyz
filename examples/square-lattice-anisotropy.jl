using Pkg
Pkg.activate(@__DIR__)

using SpinFRGLattices.SquareLattice
import PMFRG_xyz: Params, SolveFRG

J1 = 1.0
J2 = 0.5

System = getSquareLattice(6, [J1, J2])
isotropy = zeros(System.Npairs, 3)

for n = 1:System.Npairs
    isotropy[n, :] = [1.0, 0.5, 0.2]
end

Par = Params(System, N = 8, accuracy = 1e-10, temp_max = 10.0, temp_min = 1.0)

@time results = SolveFRG(Par, isotropy)
