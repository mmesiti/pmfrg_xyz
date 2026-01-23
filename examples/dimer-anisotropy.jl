using SpinFRGLattices
using JLD2
import PMFRG_xyz: Params, SolveFRG

System = SpinFRGLattices.getPolymer(2)

isotropy = zeros(System.Npairs, 3)

for n = 1:System.Npairs
    isotropy[n, :] = [1.0, 0.5, 0.2]
end


Par = Params(System, N = 8, accuracy = 1e-10, temp_max = 10.0, temp_min = 1.0)

@time sol, saved_values = SolveFRG(Par, isotropy)

save_object(
    "dimer_opt_values.jld2",
    [(saved_values.saveval[n], exp(saved_values.t[n])) for n = 1:length(saved_values.t)],
)
