using Pkg
Pkg.activate("./")

include("../src/PMFRG_xyz.jl")
using .PMFRG_xyz
include("../src/Structs.jl")

using JLD2

struct ChiAndSigma{T}
    Chi_x::Array{T,2}
    Chi_y::Array{T,2}
    Chi_z::Array{T,2}
    Sigma_x::Array{T,2}
    Sigma_y::Array{T,2}
    Sigma_z::Array{T,2}
end

function save_dynamic_chis(State, t, interator, Par)
    chi_x = PMFRG_xyz.getChi_x(State, exp(t), Par, Par.NumericalParams.N)
    chi_y = PMFRG_xyz.getChi_y(State, exp(t), Par, Par.NumericalParams.N)
    chi_z = PMFRG_xyz.getChi_z(State, exp(t), Par, Par.NumericalParams.N)
    return ObservablesDynamic(copy(chi_x), copy(chi_y), copy(chi_z))
end

function save_chi_and_sigma(State, t, integrator, Par)
    chi_x = PMFRG_xyz.getChi_x(State, exp(t), Par, Par.NumericalParams.N)
    chi_y = PMFRG_xyz.getChi_y(State, exp(t), Par, Par.NumericalParams.N)
    chi_z = PMFRG_xyz.getChi_z(State, exp(t), Par, Par.NumericalParams.N)
    return ChiAndSigma(
        copy(chi_x), copy(chi_y), copy(chi_z),
        copy(State.x[2]), copy(State.x[3]), copy(State.x[4])
    )
end

System = getSquareLattice(4, [1.0,0.5])
Par = PMFRG_xyz.Params(
    System,
    N = 6,
    temp_min = 1.0,
    temp_max = 100.0,
    accuracy = 1e-4
)

isotropy = zeros(System.Npairs, 3)
for n in 1:System.Npairs
    isotropy[n,:] .= [1.0, 0.5, 0.78]
end

let
    sol, saved_values = PMFRG_xyz.SolveFRG(
        Par, isotropy,
        SavedValues(_getFloatType(Par), ObservablesDynamic{_getFloatType(Par)}), 
        (State, t, integrator) -> save_dynamic_chis(State, t, integrator, Par)
    )
    save_object("dynamic_chi.jld2", [(saved_values.saveval[n], exp(saved_values.t[n])) for n in eachindex(saved_values.t)])
end

let
    sol, saved_values = PMFRG_xyz.SolveFRG(
        Par, isotropy,
        SavedValues(_getFloatType(Par), ChiAndSigma{_getFloatType(Par)}),
        (State, t, integrator) -> save_chi_and_sigma(State, t, integrator, Par)
    )
    save_object("chi_and_sigma.jld2", [(saved_values.saveval[n], exp(saved_values.t[n])) for n in eachindex(saved_values.t)])
end