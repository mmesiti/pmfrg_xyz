# Example script to record data for a L-flow end-to-end test

# Load dependencies from dev paths
import Pkg
Pkg.develop(path = joinpath(@__DIR__, "..", "..", "..", "..", "Recorder.jl"))
Pkg.develop(path = joinpath(@__DIR__, "..", "..", "..", "..", "SpinFRGLattices.jl"))

using JLD2
using Recorder
using SpinFRGLattices

# Include the Lambda Flow module
include("../../../src/Lpmfrg_xyz.jl")

# Create dimer system
System = SpinFRGLattices.getPolymer(2)

# Set anisotropy for all pairs
isotropy = zeros(System.Npairs, 3)
for n = 1:System.Npairs
    isotropy[n, :] = [1.0, 0.5, 0.2]
end

# Create parameters for Lambda Flow
# T = 0.5 (fixed temperature), N = 8, default lambda_min/lambda_max
Par = PMFRG_xyz.Params(
    System,
    PMFRG_xyz.LFlowNumericalParams(N = 8, T = 0.5, accuracy = 1e-10),
)

println("Recording SolveFRG for Lambda Flow regression test...")
println(
    "Parameters: N=8, T=0.5, lambda_min=$(Par.NumericalParams.lambda_min), lambda_max=$(Par.NumericalParams.lambda_max)",
)
println("System: Dimer (2 sites)")
println("Anisotropy: [1.0, 0.5, 0.2]")

# Record the SolveFRG call
@time sol, saved_values = @record PMFRG_xyz.SolveFRG(Par, isotropy)

# Create regression test data file
data_filename, script =
    Recorder.create_regression_tests(PMFRG_xyz.SolveFRG, tag = "lambda_flow")

println("\nRecording complete!")
println("Data file created: $data_filename")
println("Test script created: $script")
