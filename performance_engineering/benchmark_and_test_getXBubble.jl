using JSON
using Test
using BenchmarkTools
import ThreadPinning
import PMFRG_xyz: getXBubble!

include("benchmark_utils.jl")
include("../test/regression/dimer_anisotropy/regression_tests_dimer.jl")
include("git_utils.jl")
import .GitUtils
include("cpu_info.jl")
import .CPUInfo

function main(record::Bool, comment = "")
    println("Testing...")
    run_getXbubble_regression_tests()
    println("Benchmarking...")
    N = 6
    lattice_size = 16

    threadpinning = false

    if threadpinning
        ThreadPinning.pinthreads(:cores)
        ThreadPinning.threadinfo()
    end
    bench_result, allocations =
        benchmark_synthetic_square(N = N, lattice_size = lattice_size)

    funcnames = ["mean", "minimum", "maximum"]
    quantities = ["times", "gctimes"]
    metadata = get_metadata(comment, threadpinning, N, lattice_size)

    data = Dict(
        "benchmark_data" => Dict(
            "$(q)_$(fn)" =>
                eval(Meta.parse("$fn($(getfield(bench_result, Symbol(q))))")) for
            q in quantities for fn in funcnames
        ),
        "allocations" => allocations,
        "metadata" => metadata,
    )
    if record
        add_data_to_db(data)
    end
    println(data)
    data

end

function benchmark_synthetic_square(; N::Int = 8, lattice_size::Int = 4)
    workspace, lam = create_synthetic_workspace_square(N = N, lattice_size = lattice_size)

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

function get_metadata(comment, threadpinning, N, lattice_size)
    Dict(
        "commit" => GitUtils.get_git_commit_short(),
        "comment" => comment,
        "nthreads" => Threads.nthreads(),
        "threadpinning" => threadpinning,
        "physics" =>
            Dict("type" => "square lattice", "N" => N, "lattice_size" => lattice_size),
        "cpu_info" => CPUInfo.get_cpu_info(),
        "julia_opt_level" => Base.JLOptions().opt_level,
    )
end

function add_data_to_db(data)
    fname = joinpath(@__DIR__, "benchmark_getXBubble.db")
    if isfile(fname)
        existing_data = open(fname, "r") do f
            JSON.parse(read(f))
        end
    else
        existing_data = []
    end
    push!(existing_data, data)
    open(fname, "w") do f
        JSON.json(f, existing_data; pretty = true)
    end

end


if abspath(PROGRAM_FILE) == @__FILE__
    main(true)
end
