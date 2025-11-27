#!/usr/bin/env julia
# Benchmark script for getXBubble! function

using BenchmarkTools
using JLD2
import PMFRG_xyz: getXBubble!, addX!, addY!, ThreadLocalBuffersT
using StaticArrays
include("benchmark_utils.jl")

function main()::Int
    println("=== check addX! and addY! allocations ===\n")
    check_addXY_allocations()

    println("=== getXBubble! Benchmark ===\n")

    # # Benchmark using regression test data
    # println("1. Benchmarking with regression test data...")
    # benchmark_from_regression_data()

    # # Benchmark with synthetic data (dimer, N=8)
    # println("\n2. Benchmarking with synthetic data (dimer, N=8)...")
    # benchmark_synthetic_dimer(N=8)

    # Benchmark with synthetic data (square lattice, N=8)
    N=10
    lattice_size=8
    println("\n3. Benchmarking with synthetic data (square lattice, N=$N, lattice_size=$lattice_size)...")
    benchmark_synthetic_square(N=N, lattice_size=lattice_size)

    return 0
end

# level 1
function check_addXY_allocations()

    workspace, _ = create_synthetic_workspace_square(N=10, lattice_size=5)

    Par = workspace.Par
    (; NUnique, Npairs) = Par.System


    buffs = ThreadLocalBuffersT( zeros((21, Npairs)),
              zeros((21, Npairs)),
              zeros(21),
              zeros(3, 3, NUnique),
              zeros(3, 3, NUnique, NUnique),
              zeros(3,3),
              zeros(21),
              zeros(21),
              zeros(21),
              zeros(21))




    addX!(workspace,1,1,2,1,buffs.spropX,buffs)
    addY!(workspace,1,1,2,1,buffs.spropY,buffs)

    addXallocations = @allocations addX!(workspace,1,1,2,1,buffs.spropX,buffs)
    @assert addXallocations <= 1 "$addXallocations in addX!"

    addYallocations = @allocations addY!(workspace,1,1,2,1,buffs.spropY,buffs)
    @assert addYallocations <= 1 "$addYallocations in addY!"
end


function benchmark_from_regression_data()
    script_dir = @__DIR__
    project_root = dirname(script_dir)
    data_file = joinpath(project_root, "test/regression/dimer_anisotropy/regression_tests_dimer-PMFRG_xyz.getXBubble!.data")
    data = load_object(data_file)

    args_first = data["arguments"][1]
    workspace = args_first[1]
    lam = args_first[2]

    println("  Using first test case from regression data")
    println("  Lam = $lam")

    result = @benchmark getXBubble!($workspace, $lam)
    display(result)
    result
end

function benchmark_synthetic_dimer(; N::Int=8)
    workspace, lam = create_synthetic_workspace_dimer(N=N)

    println("  N = $N")
    println("  Lam = $lam")

    result = @benchmark getXBubble!($workspace, $lam)
    display(result)
    result
end

function benchmark_synthetic_square(; N::Int=8, lattice_size::Int=4)
    workspace, lam = create_synthetic_workspace_square(N=N, lattice_size=lattice_size)

    println("  N = $N, lattice_size = $lattice_size")
    println("  Lam = $lam")

    #getXBubble!(workspace, lam, ThreadLocalBuffers)
    #allocations = @ballocations getXBubble!($workspace, $lam, $ThreadLocalBuffers)
    #timing_results = @benchmark getXBubble!($workspace, $lam, $ThreadLocalBuffers)
    getXBubble!(workspace, lam)
    allocations = @ballocations getXBubble!($workspace, $lam)
    timing_results = @benchmark getXBubble!($workspace, $lam)
    display(allocations)
    display(timing_results)
    timing_results, allocations
end

#######

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
