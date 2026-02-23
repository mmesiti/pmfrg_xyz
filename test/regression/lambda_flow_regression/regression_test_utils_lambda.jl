function run_lambda_flow_SolveFRG_regression_tests(data_path::String)
    @testset verbose = true "Tests for PMFRG_xyz.SolveFRG (Lambda Flow)" begin
        data = load_object(data_path)
        @testset verbose = true for i = 1:length(data["return_value"])
            return_value = (data["return_value"])[i]
            arguments = (data["arguments"])[i]
            arguments_post = (data["arguments_post"])[i]
            recursive_value_test(
                return_value,
                PMFRG_xyz.SolveFRG(arguments...),
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
