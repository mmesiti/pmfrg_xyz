using Test
using SpinFRGLattices
using RecursiveArrayTools

using PMFRG_xyz
import PMFRG_xyz:
    getChi_x,
    getChi_y,
    getChi_z,
    getChi_x_refactored,
    getChi_y_refactored,
    getChi_z_refactored,
    _getFloatType

function make_random_tflow_state(Par)
    NUnique = Par.System.NUnique
    N = Par.NumericalParams.N
    VDims = PMFRG_xyz.getVDims(Par)
    FT = _getFloatType(Par)
    return ArrayPartition(
        randn(FT, NUnique),
        randn(FT, NUnique, N),
        randn(FT, NUnique, N),
        randn(FT, NUnique, N),
        randn(FT, VDims...),
    )
end

function test_chi_refactored_vs_original_tflow()
    System = SpinFRGLattices.getPolymer(2)
    Par = PMFRG_xyz.Params(System, N = 4, accuracy = 1e-6, temp_min = 0.5, temp_max = 2.0)
    T = 1.0
    @testset "getChi_z_refactored vs getChi_z (T-flow)" begin
        for _ = 1:3
            State = make_random_tflow_state(Par)
            @test getChi_z_refactored(State, T, Par) ≈ getChi_z(State, T, Par)
        end
    end
    @testset "getChi_x_refactored vs getChi_x (T-flow)" begin
        for _ = 1:3
            State = make_random_tflow_state(Par)
            @test getChi_x_refactored(State, T, Par) ≈ getChi_x(State, T, Par)
        end
    end
    @testset "getChi_y_refactored vs getChi_y (T-flow)" begin
        for _ = 1:3
            State = make_random_tflow_state(Par)
            @test getChi_y_refactored(State, T, Par) ≈ getChi_y(State, T, Par)
        end
    end
end

test_chi_refactored_vs_original_tflow()
