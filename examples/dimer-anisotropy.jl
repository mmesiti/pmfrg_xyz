using Pkg
Pkg.activate(@__DIR__)

using SpinFRGLattices
import PMFRG_xyz: Params, SolveFRG

System = SpinFRGLattices.getPolymer(2)
par = Params(System, N = 8, accuracy = 1e-10, temp_max = 10.0, temp_min = 1.0)
isotropy = zeros(System.Npairs, 3)

for n = 1:System.Npairs
    isotropy[n, :] = [1.0, 0.5, 0.2]
end


@time results = SolveFRG(par, isotropy)
