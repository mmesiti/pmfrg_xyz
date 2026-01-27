# PMFRG_xyz-specific regression test utilities
# This file contains test functions specific to the T Flow (PMFRG_xyz module)

using PMFRG_xyz

# Include generic utilities that work for both T Flow and Lambda Flow
include("regression_test_utils_generic.jl")

# Parameterized regression test functions for PMFRG_xyz (T Flow)

function run_getXbubble_regression_tests(data_path::String)
    @testset verbose = true "Tests for PMFRG_xyz.getXBubble!" begin
        data = load_object(data_path)
        @testset verbose = false for i = 1:length(data["return_value"])
            return_value = (data["return_value"])[i]
            arguments = (data["arguments"])[i]
            arguments_post = (data["arguments_post"])[i]
            recursive_value_test(
                return_value,
                (PMFRG_xyz.getXBubble!)(arguments...),
                "return values - case $i",
                true,
            )
            @testset "arguments" begin
                for i in eachindex(arguments)
                    recursive_value_test(arguments[i], arguments_post[i], "idx = $i", true)
                end
            end
        end

        # Test Float32 precision parameter
        @testset verbose = true "Float32 precision test" begin
            arguments = (data["arguments"])[1]
            Workspace = arguments[1]
            T = arguments[2]

            X_original = copy(Workspace.X)

            PMFRG_xyz.setZero!(Workspace.X)
            @test_nowarn PMFRG_xyz.getXBubble!(Workspace, T)
            X_float64 = copy(Workspace.X)

            PMFRG_xyz.setZero!(Workspace.X)
            @test_nowarn PMFRG_xyz.getXBubble!(Workspace, T; ComputeType = Float64)
            X_float64_explicit = copy(Workspace.X)

            PMFRG_xyz.setZero!(Workspace.X)
            @test_nowarn PMFRG_xyz.getXBubble!(Workspace, T; ComputeType = Float32)
            X_float32 = copy(Workspace.X)

            @test X_float64 ≈ X_float64_explicit
            @test X_float64 ≈ X_float32 rtol = 1e-6

            Workspace.X .= X_original
        end
    end
    nothing
end

function run_SolveFRG_regression_tests(data_path::String)
    @testset verbose = true "Tests for PMFRG_xyz.SolveFRG" begin
        data = load_object(data_path)
        @testset verbose = true for i = 1:length(data["return_value"])
            return_value = (data["return_value"])[i]
            arguments = (data["arguments"])[i]
            arguments_post = (data["arguments_post"])[i]
            recursive_value_test(
                return_value,
                (PMFRG_xyz.SolveFRG)(arguments...),
                "return values - case $i",
                true,
            )
            @testset verbose = true for i in eachindex(arguments)
                recursive_value_test(
                    arguments[i],
                    arguments_post[i],
                    "arguments - case $i",
                    false,
                )
            end
        end
    end
end
