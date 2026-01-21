using Test
using PMFRG_xyz

include("regression/dimer_anisotropy/regression_tests_dimer.jl")
include("performance/allocations.jl")
include("unit/vertex_functions.jl")

function run_regression_tests()
    @testset verbose = true "Regression Tests for PMFRG_xyz, dimer + anisotropy" begin
        run_getXbubble_regression_tests()
        run_SolveFRG_regression_tests()
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
        @testset "Vertex function equivalence tests" begin
            test_fillvbuffer_v_equivalence()
        end
    end
end

@time @testset verbose = true "PMFRG_xyz tests" begin
    run_unit_tests()
    run_regression_tests()
    run_allocation_tests()
end
