# Performance Engineering for pmfrg_xyz

This directory contains tools for profiling and benchmarking pmfrg_xyz.

## Installed Packages

- **PProf.jl**: Statistical profiling with pprof integration
- **BenchmarkTools.jl**: Accurate benchmarking framework
- **ProfileSVG.jl**: SVG-based profile visualization
- **FlameGraphs.jl**: Flame graph visualization
- **LIKWID.jl**: Hardware performance counter access
- **ThreadPinning.jl**: Thread affinity control

## Profiling

### PProf Statistical Profiling

Profile examples with PProf (`dimer` or `square_lattice`):
```bash
julia pmfrg_xyz/performance_engineering/profile_pprof.jl [example_name]
```

Examples:
```bash
# Profile dimer example (default)
julia pmfrg_xyz/performance_engineering/profile_pprof.jl
julia pmfrg_xyz/performance_engineering/profile_pprof.jl dimer

# Profile square lattice example
julia pmfrg_xyz/performance_engineering/profile_pprof.jl square_lattice
```

This will:
- Run the calculation once for compilation
- Profile the second run with detailed call stack information
- Display PProf version and executable path
- Save profile data to `profile_data/profile_<example>_<commit>_pprof<version>.pb.gz`
- Open an interactive pprof viewer in your browser

### LIKWID Hardware Counter Profiling

Profile specific functions with hardware performance counters:
```bash
# Apply instrumentation patch
git apply performance_engineering/likwid_markers.patch
julia --project=. -e 'using Pkg; Pkg.add("LIKWID")'

# Run profiling (single-threaded, pinned to core 0)
cd performance_engineering
./run_likwid_profile_getXBubble.sh FLOPS_DP  

# Revert instrumentation
cd ..
git restore ./src/PMFRG_xyz.jl Manifest.toml Project.toml
```
At the moment only performance data for the the components of the getXBubble function is collected.

**Instrumented functions:**
- `launchPMFRG!` - main FRG solver
- `getDeriv!` - derivative computation orchestrator
- `getXBubble!` - bubble diagram computation
- `addX!` - X-channel contributions
- `addY!` - Y-channel contributions

List available groups: `likwid-perfctr -a`

## Benchmarking

### getXBubble! Benchmark and Test

Combined benchmarking and regression testing for the `getXBubble!` function:

```bash
julia --project=. performance_engineering/benchmark_and_test_getXBubble.jl
```

This script:
1. **Runs regression tests** to ensure correctness
2. **Benchmarks `getXBubble!`** with square lattice (`N=6`, `lattice_size=16`)
3. **Records results** to `benchmark_getXBubble.db` (JSON format)

**Tracked data:**
- Timing statistics: mean, minimum, maximum (times and GC times)
- Allocation counts
- Git commit hash
- CPU model and frequency
- Thread count and pinning status
- Physics parameters (`N`, `lattice_size`)

**Customization:** Edit the script to modify:
- `N` (frequency discretization, default: 6)
- `lattice_size` (system size, default: 16)
- `threadpinning` (enable/disable, default: false)

## Utilities

### Helper Modules

- **`benchmark_utils.jl`**: Synthetic workspace creation for benchmarks
  - `create_synthetic_workspace_dimer(N)`
  - `create_synthetic_workspace_square(N, lattice_size)`
  - `check_addXY_allocations()` - verify zero-allocation property

- **`example_configs.jl`**: Example system configurations
  - Dimer and square lattice parameter sets

- **`git_utils.jl`**: Git information utilities
  - `get_git_commit_short()` - current commit hash
  - `check_git_status()` - verify clean working tree

- **`cpu_info.jl`**: CPU information extraction
  - Reads `/proc/cpuinfo` and SLURM environment
  - Reports model name, frequency, core count

- **`pprof_version.jl`**: PProf version utilities (used by `profile_pprof.jl`)
  - `get_pprof_version()` - extract pprof version string
  - `get_pprof_path()` - get pprof executable path

### SLURM Integration

**`slurm-profiling.sh`**: Stub script for SLURM cluster profiling jobs
