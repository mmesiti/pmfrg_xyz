using JLD2
using SpinFRGLattices
using Test
# Include the Lambda Flow module
include("../../../src/PMFRG_xyz.jl")


using .PMFRG_xyz
using RecursiveArrayTools

"""
Function that can be used to generate a random state object for testing.
"""
function make_random_state_and_par()
    System = SpinFRGLattices.getPolymer(2)
    Par = PMFRG_xyz.Params(
        System,
        PMFRG_xyz.LFlowNumericalParams(
            N = 4,
            accuracy = 1e-6,
            T = 0.4,
            lambda_min = 0.5,
            lambda_max = 2.0,
        ),
    )


    NUnique = Par.System.NUnique
    N = Par.NumericalParams.N
    VDims = PMFRG_xyz.getVDims(Par)
    FT = PMFRG_xyz._getFloatType(Par)
    return ArrayPartition(
        randn(FT, NUnique),
        randn(FT, NUnique, N),
        randn(FT, NUnique, N),
        randn(FT, NUnique, N),
        randn(FT, VDims...),
    ),
    Par
end


function get_data()
    (; State, FlowParam, Par) =
        load_object(joinpath(@__DIR__, "get_Chi_dimer_random_input.jld2"))
    (; chi_x, chi_y, chi_z) = load_object(joinpath(@__DIR__, "get_Chi_dimer_output.jld2"))
    return (; State, FlowParam, Par, chi_x, chi_y, chi_z)
end


function test_chi_Lflow()
    (; State, FlowParam, Par, chi_x, chi_y, chi_z) = get_data()

    @testset "L-Flow" begin
        @test PMFRG_xyz.getChi_z(State, FlowParam, Par) ≈ chi_z
        @test PMFRG_xyz.getChi_x(State, FlowParam, Par) ≈ chi_x
        @test PMFRG_xyz.getChi_y(State, FlowParam, Par) ≈ chi_y
    end
end
