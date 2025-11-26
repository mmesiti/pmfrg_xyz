using JSON


include("benchmark_getXBubble.jl")
include("../test/regression/dimer_anisotropy/regression_tests_dimer.jl")
include("git_utils.jl")
import .GitUtils
include("cpu_info.jl")
import .CPUInfo
import ThreadPinning

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
        JSON.json(f, existing_data; pretty=true)
    end

end

function getXBubble_test_and_benchmark(record::Bool, comment="")
    println("Testing...")
    run_getXbubble_regression_tests()
    println("Benchmarking...")
    N = 10
    lattice_size = 5

    threadpinning = false 

    if threadpinning 
        ThreadPinning.pinthreads(:cores)
        ThreadPinning.threadinfo()
    end
    bench_result, allocations = benchmark_synthetic_square(N=N, lattice_size=lattice_size)

    funcnames = ["mean", "minimum", "maximum"]
    quantities = ["times", "gctimes"]


    data = Dict("commit" => GitUtils.get_git_commit_short(),
        "comment" => comment,
        "nthreads" => Threads.nthreads(),
        "threadpinning" => threadpinning,
        "benchmark_data" => Dict("$(q)_$(fn)" => eval(Meta.parse("$fn($(getfield(bench_result, Symbol(q))))"))
                                 for q in quantities
                                 for fn in funcnames),
        "allocations" => allocations,
        "physics" => Dict(
            "type" => "square lattice",
            "N" => N,
            "lattice_size" => lattice_size,
        ),
        "cpu_info" => CPUInfo.get_cpu_info())
    if record
        add_data_to_db(data)
    end
    return data
end



if abspath(PROGRAM_FILE) == @__FILE__
    println(getXBubble_test_and_benchmark(true))
end
