using Test
using PMFRG_xyz

include("regression/dimer_anisotropy/regression_tests_dimer.jl")
include("performance/allocations.jl")
include("unit/Xtype.jl")

function run_regression_tests()
    @testset verbose = true "Regression Tests for PMFRG_xyz, dimer anisotropy" begin
        run_getXbubble_regression_tests()
        run_SolveFRG_regression_tests()
    end
end

function run_allocation_tests()
    @testset "Checking allocations" begin
        check_addXY_allocations()
    end
end

function run_unit_tests()
    @testset verbose = true "Unit Tests" begin
        @testset "XType tests" begin
            # Test with small configuration
            mapping_small = DefaultXIndexMapping(21, 3, 4)
            run_all_tests(mapping_small)

            # Test with realistic configuration
            mapping_realistic = DefaultXIndexMapping(42, 10, 24)
            run_all_tests(mapping_realistic)
        end
    end
end

@time @testset verbose = true "PMFRG_xyz tests" begin
    run_unit_tests()
    run_regression_tests()
    run_allocation_tests()
end
