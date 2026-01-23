import Pkg
Pkg.activate(".")

using CairoMakie

### Dimer T-Flow ###

include("Tpmfrg_xyz.jl")
using .Tpmfrg_xyz
using SpinFRGLattices
using OrdinaryDiffEq
using JLD2

System = getPolymer(2)
Par = Tpmfrg_xyz.Params(System, N = 8, accuracy = 1e-5, temp_min = 0.1, temp_max = 100.0)

isotropy = zeros(System.Npairs, 3)
for n = 1:System.Npairs
    isotropy[n, :] = [0.4, 0.0, 1.0]
end

isotropy

sol, saved_values = Tpmfrg_xyz.SolveFRG(Par, isotropy, method = DP5(), save_steps = true)
save_object(
    "TFlow_dimer.jld2",
    [(saved_values.saveval[n], exp(saved_values.t[n])) for n = 1:length(saved_values.t)],
)

####################


function ChiExact12(T, w, jx, jy, jz)
    b = 1 / T
    Jp = (jy + jx) / 2.0
    Jm = (jx - jy) / 2.0
    Jyz = (jz + jy) / 2.0
    Jxz = (jz + jx) / 2.0

    v = 0
    if (Jm == 0)
        v = exp(-10)
    end

    Z = 2.0 * (1.0 + exp(-b * Jp) + exp(-b * Jxz) + exp(-b * Jyz))
    A = -(1.0 - exp(-b * Jp)) / Jp
    B = -(1.0 - exp(b * Jm)) * exp(-b * Jxz) / Jm

    if (isapprox(Jm, 0.0, atol = 1e-5))
        B = (b + b^2 * Jm) * exp(-b * Jxz)
    end

    result = (A + B) / Z

    return result
end

function ChiExact11(T, w, jx, jy, jz)
    b = 1 / T
    Jp = (jy + jx) / 2.0
    Jm = (jx - jy) / 2.0
    Jyz = (jz + jy) / 2.0
    Jxz = (jz + jx) / 2.0

    v = 0
    if (Jm == 0)
        v = exp(-10)
    end

    Z = 2.0 * (1.0 + exp(-b * Jp) + exp(-b * Jxz) + exp(-b * Jyz))
    A = (1.0 - exp(-b * Jp)) / Jp
    B = -(1.0 - exp(b * Jm)) * exp(-b * Jxz) / Jm

    if (isapprox(Jm, 0.0, atol = 1e-5))
        B = (b + b^2 * Jm) * exp(-b * Jxz)
    end

    result = (A + B) / Z

    return result
end

Chi_z(T, jx, jy, jz) = ChiExact12(T, 0.0, jx, jy, jz)
Chi_x(T, jx, jy, jz) = ChiExact12(T, 0.0, jz, jy, jx)
Chi_y(T, jx, jy, jz) = ChiExact12(T, 0.0, jx, jz, jy)
ChiLocal_z(T, jx, jy, jz) = ChiExact11(T, 0.0, jx, jy, jz)
ChiLocal_x(T, jx, jy, jz) = ChiExact11(T, 0.0, jz, jy, jx)
ChiLocal_y(T, jx, jy, jz) = ChiExact11(T, 0.0, jx, jz, jy)

function ChiHom(T, w, j)
    b = 1 / T
    return (-(exp(b) - 1 - b) / (2 * (exp(b) + 3)) * (1 / (1 + w^2)))
end

chiVals_x = []
chiVals_y = []
chiVals_z = []

using JLD2
data = load_object("TFlow_dimer.jld2")

length(data)

for n in eachindex(data)
    append!(chiVals_x, data[n][1].Chi_x[2])
    append!(chiVals_y, data[n][1].Chi_y[2])
    append!(chiVals_z, data[n][1].Chi_z[2])
end

isotropy = [0.4, 0.0, 1.0]

T = [data[n][end] for n in eachindex(data)]

fig = Figure()
ax = Axis(fig[1, 1], ylabel = L"χ", xlabel = L"T", title = "T-Flow J=[0.4,0,1]")
n = 400
lines!(
    ax,
    T[n:end],
    Chi_x.(T[n:end], isotropy[1], isotropy[2], isotropy[3]),
    linewidth = 2,
    label = L"χ_x",
)
lines!(
    ax,
    T[n:end],
    Chi_y.(T[n:end], isotropy[1], isotropy[2], isotropy[3]),
    linewidth = 2,
    label = L"χ_y",
)
lines!(
    ax,
    T[n:end],
    Chi_z.(T[n:end], isotropy[1], isotropy[2], isotropy[3]),
    linewidth = 2,
    label = L"χ_z",
)
scatter!(ax, T[n:end], chiVals_x[n:end], marker = :diamond, markersize = 5)
scatter!(ax, T[n:end], chiVals_y[n:end], marker = :diamond, markersize = 5)
scatter!(ax, T[n:end], chiVals_z[n:end], marker = :diamond, markersize = 5)
# lines!(ax, x, ChiHom.(x, 0, 1.0))
axislegend(ax, position = :rt)
display("image/png", fig)
