#!/usr/bin/env julia
# Usage: julia profile.jl [example_name]

using Pkg
Pkg.activate(@__DIR__)

using PProf, Profile, PMFRG_xyz
import ThreadPinning

include("example_configs.jl")
using .ExampleSetups: example_setups

include("git_utils.jl")
using .GitUtils: get_git_commit_short, check_git_status

include("pprof_version.jl")
using .PProfVersionUtils: get_pprof_version, get_pprof_path

# level 0
function main()
    setup_output_dir()
    example_name = parse_args()
    check_git_status()

    par, isotropy = setup_example(example_name)
    run_compilation(par, isotropy)
    run_profiling(par, isotropy)
    save_and_display_results(example_name)
end

# level 1
function setup_output_dir()
    mkpath(joinpath(@__DIR__, "profile_data"))
end

function parse_args()
    length(ARGS) >= 1 ? ARGS[1] : "dimer"
end

function setup_example(example_name)
    println("Setting up $example_name example...")
    example_setups[example_name]()
end

function run_compilation(par, isotropy)
    println("Running initial execution (for compilation)...")
    SolveFRG(par, isotropy)
end

function run_profiling(par, isotropy)
    println("\nStarting profiling run...")
    ThreadPinning.pinthreads(:cores)
    ThreadPinning.threadinfo()
    Profile.clear()
    @profile SolveFRG(par, isotropy)
end

function save_and_display_results(example_name)
    git_commit = get_git_commit_short()
    pprof_version = get_pprof_version()
    pprof_path = get_pprof_path()
    profile_file = get_profile_filename(example_name, git_commit, pprof_version)

    println("\nSaving profile data to: $profile_file")
    pprof(out = profile_file)

    println("\nGenerating interactive profile visualization...")
    pprof()

    print_summary(profile_file, git_commit, pprof_version, pprof_path, example_name)
end

# level 2
function get_profile_filename(example_name, git_commit, pprof_version)
    slurm_job_id=get(ENV,"SLURM_JOB_ID","noslurm")
    slurm_cpu_freq_req=get(ENV,"SLURM_CPU_FREQ_REQ","nofreq")
    joinpath(
        @__DIR__,
        "profile_data",
        "profile_$(example_name)_" *
        "$(git_commit)_" * 
        "pprof$(pprof_version)_" * 
        "$(slur_job_id)_" *
        "$(slurm_cpu_freq_req).pb.gz",
    )
end

function print_summary(profile_file, git_commit, pprof_version, pprof_path, example_name)
    println("\nProfile data saved to: $profile_file")
    println("Git commit: $git_commit")
    println("PProf version: $pprof_version")
    println("PProf executable: $pprof_path")
    println("Example: $example_name")
    println("\nTo view saved profile later, run:")
    println("  julia -e 'using PProf; pprof(\"$profile_file\")'")
    println("Or use the pprof command-line tool:")
    println("  $pprof_path -http=: \"$profile_file\"")
end

#######

main()
