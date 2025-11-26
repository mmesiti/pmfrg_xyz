#!/usr/bin/env julia
# Benchmark script for getXBubble! function

using BenchmarkTools
using JLD2
import PMFRG_xyz: getXBubble!, addX!, addY!
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

    workspace, _ = create_synthetic_workspace_dimer(N=10)

    Par = workspace.Par
    (; NUnique) = Par.System


    buffs = (V12=zeros((21, maximum(Par.System.siteSum.ki))),
        V34=zeros((21, maximum(Par.System.siteSum.kj))),
        X_sum=(@MVector zeros(21)),
        spropX=(@MArray zeros(3, 3, NUnique)),
        spropY=zeros(3, 3, NUnique, NUnique))

    (;spropX,spropY) = buffs

    spropX = (@MArray zeros(3,3,NUnique))

    addX!(workspace,1,1,2,1,spropX,buffs)
    addY!(workspace,1,1,2,1,spropY)

    addXallocated = @allocated addX!(workspace,1,1,2,1,spropX,buffs)
    addXallocations = @allocations addX!(workspace,1,1,2,1,spropX,buffs)
    addYallocated = @allocated addY!(workspace,1,1,2,1,spropY)
    addYallocations = @allocations addY!(workspace,1,1,2,1,spropY)

    println("addXAllocations: $addXallocations ($addXallocated)")
    println("addYAllocations: $addYallocations ($addYallocated)")
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
