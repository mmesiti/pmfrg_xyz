#!/usr/bin/env julia
# Benchmark script for getXBubble! function

using BenchmarkTools
using JLD2
using SpinFRGLattices
using SpinFRGLattices.SquareLattice
import PMFRG_xyz: Params, getXBubble!, AllocateSetup, InitializeState, OneLoopWorkspace

function main()::Int
    println("=== getXBubble! Benchmark ===\n")

    # Benchmark using regression test data
    println("1. Benchmarking with regression test data...")
    benchmark_from_regression_data()

    # Benchmark with synthetic data (dimer, N=8)
    println("\n2. Benchmarking with synthetic data (dimer, N=8)...")
    benchmark_synthetic_dimer(N=8)

    # Benchmark with synthetic data (square lattice, N=8)
    N=10
    lattice_size=8
    println("\n3. Benchmarking with synthetic data (square lattice, N=$N, lattice_size=$lattice_size)...")
    benchmark_synthetic_square(N=N, lattice_size=lattice_size)

    return 0
end

# level 1
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

    allocations = @ballocations getXBubble!($workspace, $lam)
    timing_results = @benchmark getXBubble!($workspace, $lam)
    display(timing_results)
    timing_results, allocations
end

# level 2
function create_synthetic_workspace_dimer(; N::Int=8)
    System = SpinFRGLattices.getPolymer(2)
    par = Params(System, N=N, temp_max=10.0, temp_min=1.0)
    isotropy = zeros(System.Npairs, 3)

    for n in 1:System.Npairs
        isotropy[n, :] = [1.0, 0.5, 0.2]
    end

    State = InitializeState(par, isotropy)
    setup = AllocateSetup(par)

    X = setup[1]
    Par = setup[end]
    Deriv = copy(State)
    fill_with_zeros!(Deriv)

    workspace = OneLoopWorkspace(State, Deriv, X, Par)
    lam = 5.0

    return workspace, lam
end

function create_synthetic_workspace_square(; N::Int=8, lattice_size::Int=4)
    J1 = 1.0
    J2 = 0.5

    System = getSquareLattice(lattice_size, [J1, J2])
    par = Params(System, N=N, temp_max=10.0, temp_min=1.0)
    isotropy = zeros(System.Npairs, 3)

    for n in 1:System.Npairs
        isotropy[n, :] = [1.0, 0.5, 0.2]
    end

    State = InitializeState(par, isotropy)
    setup = AllocateSetup(par)

    X = setup[1]
    Par = setup[end]
    Deriv = copy(State)
    fill_with_zeros!(Deriv)

    workspace = OneLoopWorkspace(State, Deriv, X, Par)
    lam = 5.0

    return workspace, lam
end

# level 3
function fill_with_zeros!(state)
    if hasfield(typeof(state), :f_int)
        fill!(state.f_int, 0.0)
    end
    if hasfield(typeof(state), :iSigma)
        if hasfield(typeof(state.iSigma), :x)
            fill!(state.iSigma.x, 0.0)
        end
        if hasfield(typeof(state.iSigma), :y)
            fill!(state.iSigma.y, 0.0)
        end
        if hasfield(typeof(state.iSigma), :z)
            fill!(state.iSigma.z, 0.0)
        end
    end
    if hasfield(typeof(state), :Gamma)
        fill!(state.Gamma, 0.0)
    end
end

#######

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
