import PMFRG_xyz: Params, SolveFRG, TFlow
using JLD2
using SpinFRGLattices.SquareLattice

J1 = 1.0
J2 = 0.5

System = getSquareLattice(6, [J1, J2])
isotropy = zeros(System.Npairs, 3)

for n = 1:System.Npairs
    isotropy[n, :] = [1.0, 0.5, 0.2]
end

Par = Params(System, TFlow(), N = 8, accuracy = 1e-10, temp_max = 10.0, temp_min = 1.0)

@time sol, saved_values = SolveFRG(Par, isotropy)

save_object(
    "square_lattice_values.jld2",
    [(saved_values.saveval[n], exp(saved_values.t[n])) for n = 1:length(saved_values.t)],
)
