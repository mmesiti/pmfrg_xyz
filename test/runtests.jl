using Test
using PMFRG_xyz

include("regression/dimer_anisotropy/regression_tests_dimer.jl")
include("performance/allocations.jl")

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

@time @testset verbose=true "PMFRG_xyz tests" begin
    run_regression_tests()
    run_allocation_tests()
end

