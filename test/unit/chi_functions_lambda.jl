using Test
using RecursiveArrayTools
using SpinFRGLattices

include("../../src/Lpmfrg_xyz.jl")
using .pmfrg_xyz

function make_random_lflow_state(Par)
    NUnique = Par.System.NUnique
    N = Par.NumericalParams.N
    VDims = pmfrg_xyz.getVDims(Par)
    FT = pmfrg_xyz._getFloatType(Par)
    return ArrayPartition(
        randn(FT, NUnique),
        randn(FT, NUnique, N),
        randn(FT, NUnique, N),
        randn(FT, NUnique, N),
        randn(FT, VDims...),
    )
end

function test_chi_refactored_vs_original_lambda_flow()
    System = SpinFRGLattices.getPolymer(2)
    Par = pmfrg_xyz.Params(System, N = 4, T = 0.5, accuracy = 1e-6)
    Lam = 1.0
    @testset "getChi_z_refactored vs getChi_z (lambda-flow)" begin
        for _ = 1:3
            State = make_random_lflow_state(Par)
            @test pmfrg_xyz.getChi_z_refactored(State, Lam, Par) ≈
                  pmfrg_xyz.getChi_z(State, Lam, Par)
        end
    end
    @testset "getChi_x_refactored vs getChi_x (lambda-flow)" begin
        for _ = 1:3
            State = make_random_lflow_state(Par)
            @test pmfrg_xyz.getChi_x_refactored(State, Lam, Par) ≈
                  pmfrg_xyz.getChi_x(State, Lam, Par)
        end
    end
    @testset "getChi_y_refactored vs getChi_y (lambda-flow)" begin
        for _ = 1:3
            State = make_random_lflow_state(Par)
            @test pmfrg_xyz.getChi_y_refactored(State, Lam, Par) ≈
                  pmfrg_xyz.getChi_y(State, Lam, Par)
        end
    end
end

test_chi_refactored_vs_original_lambda_flow()
