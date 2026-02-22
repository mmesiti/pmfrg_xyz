# Lambda Flow (pmfrg_xyz) regression test runner

using SpinFRGLattices

# Include the Lambda Flow module
include("../../../src/Lpmfrg_xyz.jl")

# include chi regression tests
include("./chi-dimer.jl")

# Include generic regression test utilities
include("../regression_test_utils_generic.jl")

# Run Lambda Flow SolveFRG regression tests
function run_lambda_flow_SolveFRG_regression_tests(data_path::String)
    @testset verbose = true "Tests for pmfrg_xyz.SolveFRG (Lambda Flow)" begin
        data = load_object(data_path)
        @testset verbose = true for i = 1:length(data["return_value"])
            return_value = (data["return_value"])[i]
            arguments = (data["arguments"])[i]
            arguments_post = (data["arguments_post"])[i]
            recursive_value_test(
                return_value,
                pmfrg_xyz.SolveFRG(arguments...),
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
        @testset verbose = true "chi" begin
            test_chi_Lflow()
        end
    end
end

# Run the tests
data_file = joinpath(@__DIR__, "regression_tests_lambda_flow.data")
println("Loading regression test data from: $data_file")
run_lambda_flow_SolveFRG_regression_tests(data_file)

println("\nRegression tests complete!")
