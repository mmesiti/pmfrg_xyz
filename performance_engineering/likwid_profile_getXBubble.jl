#!/usr/bin/env julia
# LIKWID profiling script for PMFRG_xyz
# Usage: likwid-perfctr -C 0 -g GROUP -m julia --project=. -t 1 likwid_profile.jl [GROUP]

using LIKWID
import ThreadPinning
using Pkg
using ArgParse
Pkg.instantiate()

include("benchmark_utils.jl")
import PMFRG_xyz: getXBubble!, getDeriv!, SolveFRG

function main(precision)
    N = 10
    lattice_size = 16

    setup_threads()
    run_profiling(N, lattice_size, precision)
    return 0
end

# level 1

function setup_threads()
    ThreadPinning.pinthreads(:cores)
    println("Thread pinning info:")
    ThreadPinning.threadinfo()
end

function run_profiling(N::Int, lattice_size::Int, precision::Int)
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

    T = Dict(64 => Float64, 32 => Float32)
    @marker "getXBubble!" getXBubble!(workspace, lam, ComputeType = T[precision])

    Marker.close()
    println("Done.")
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--precision"
        help = "64 -> use FP64, 32 -> use FP32"
        arg_type = Int
        default = 64
    end

    return parse_args(s)
end


if abspath(PROGRAM_FILE) == @__FILE__
    args = parse_commandline()
    exit(main(args["precision"]))
end
