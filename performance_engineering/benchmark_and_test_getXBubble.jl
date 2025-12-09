using JSON
using Test
using BenchmarkTools
using ArgParse
import ThreadPinning
import PMFRG_xyz: getXBubble!

include("benchmark_utils.jl")
include("../test/regression/dimer_anisotropy/regression_tests_dimer.jl")
include("git_utils.jl")
import .GitUtils
include("cpu_info.jl")
import .CPUInfo

function main(;
    record::Bool,
    commit::Union{String,Nothing} = nothing,
    dbfile::String,
    comment::String,
)
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
    metadata = get_metadata(comment, threadpinning, N, lattice_size, commit)

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
        add_data_to_db(data, dbfile)
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

function get_metadata(
    comment,
    threadpinning,
    N,
    lattice_size,
    commit::Union{String,Nothing} = nothing,
)
    commit_hash = isnothing(commit) ? GitUtils.get_git_commit_short() : commit
    Dict(
        "commit" => commit_hash,
        "comment" => comment,
        "nthreads" => Threads.nthreads(),
        "threadpinning" => threadpinning,
        "physics" =>
            Dict("type" => "square lattice", "N" => N, "lattice_size" => lattice_size),
        "cpu_info" => CPUInfo.get_cpu_info(),
        "julia_opt_level" => Base.JLOptions().opt_level,
    )
end

function add_data_to_db(data, dbfilename)
    fname = dbfilename
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


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--no-save"
        help = "do not save results to database"
        action = :store_true
        "--commit"
        help = "git commit hash (default: auto-detect current commit)"
        arg_type = String
        default = nothing
        "--comment"
        help = "comment to add to the recording"
        arg_type = String
        default = ""
        "--dbfile"
        help = "file where to save the benchmark results and metadata"
        arg_type = String
        default = joinpath(@__DIR__, "benchmark_getXBubble.db")
    end

    return parse_args(s)
end

if abspath(PROGRAM_FILE) == @__FILE__
    args = parse_commandline()
    save_to_db = !args["no-save"]
    commit = args["commit"]
    comment = args["comment"]
    dbfilename = args["dbfile"]
    main(; record = save_to_db, commit = commit, dbfile = dbfilename, comment = comment)
end
