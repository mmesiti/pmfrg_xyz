using Test
using PMFRG_xyz

include("regression/regression_test_utils.jl")
include("regression/chi/dimer.jl")
include("regression/lambda_flow_regression/regression_test_utils_lambda.jl")
include("regression/lambda_flow_regression/chi-dimer.jl")
include("performance/allocations.jl")
include("unit/vertex_functions.jl")

const DIMER_DIR = joinpath(@__DIR__, "regression", "dimer_anisotropy")
const SQUARE_LATTICE_DIR = joinpath(@__DIR__, "regression", "square_lattice_anisotropy")
const LAMBDA_FLOW_DIR = joinpath(@__DIR__, "regression", "lambda_flow_regression")

function run_regression_tests()
    @testset verbose = true "T-flow" begin
        @testset verbose = true "getChi" begin
            test_chi_Tflow()
        end
        @testset verbose = true "Regression Tests for PMFRG_xyz, dimer + anisotropy" begin
            run_getXbubble_regression_tests(
                joinpath(DIMER_DIR, "regression_tests_dimer-PMFRG_xyz.getXBubble!.data"),
            )
            run_SolveFRG_regression_tests(
                joinpath(DIMER_DIR, "regression_tests_dimer-PMFRG_xyz.SolveFRG.data"),
            )
        end
        @testset verbose = true "Regression Tests for PMFRG_xyz, square lattice + anisotropy" begin
            run_getXbubble_regression_tests(
                joinpath(
                    SQUARE_LATTICE_DIR,
                    "regression_tests_square_lattice-PMFRG_xyz.getXBubble!.data",
                ),
            )
            run_SolveFRG_regression_tests(
                joinpath(
                    SQUARE_LATTICE_DIR,
                    "regression_tests_square_lattice-PMFRG_xyz.SolveFRG.data",
                ),
            )
        end

    end

    @testset verbose = true "L-flow" begin
        @testset verbose = true "getChi" begin
            test_chi_Lflow()
        end
        @testset verbose = true "Regression Tests for PMFRG_xyz, Lambda Flow, dimer" begin
            run_lambda_flow_SolveFRG_regression_tests(
                joinpath(LAMBDA_FLOW_DIR, "regression_tests_lambda_flow.data"),
            )
        end


    end
end

function run_allocation_tests()
    @testset "Checking allocations" for T in (Float64, Float32)
        addXallocations, addYallocations = check_addXY_allocations(T)
        @test addXallocations <= 2
        @test addYallocations <= 2
    end
end

function run_unit_tests()
    @testset verbose = true "Unit Tests" begin
        @testset "Equivalence between FillVBuffer! and V_" begin
            test_fillvbuffer_v_equivalence()
        end
    end
end

@time @testset verbose = true "PMFRG_xyz tests" begin
    run_unit_tests()
    run_regression_tests()
    run_allocation_tests()
end
