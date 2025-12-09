#!/usr/bin/env julia
# LIKWID profiling script for PMFRG_xyz
# Usage: likwid-perfctr -C 0 -g GROUP -m julia --project=. -t 1 likwid_profile.jl [GROUP]

using LIKWID
import ThreadPinning
using Pkg
Pkg.instantiate()

include("benchmark_utils.jl")
import PMFRG_xyz: getXBubble!, getDeriv!, SolveFRG

function main()
    N = 10
    lattice_size = 16

    setup_threads()
    run_profiling(N, lattice_size)
    return 0
end

# level 1

function setup_threads()
    ThreadPinning.pinthreads(:cores)
    println("Thread pinning info:")
    ThreadPinning.threadinfo()
end

function run_profiling(N::Int, lattice_size::Int)
    println("\n=== LIKWID Profiling ===")
    println("N=$N, lattice_size=$lattice_size")
    println("Threads: $(Threads.nthreads())")

    workspace, lam = create_synthetic_workspace_square(N = N, lattice_size = lattice_size)

    # Warmup run (outside markers)
    println("\nWarmup...")
    getXBubble!(workspace, lam)

    # Profiled run
    println("Profiling getXBubble!...")
    Marker.init()

    @marker "getXBubble!" getXBubble!(workspace, lam)

    Marker.close()
    println("Done.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
