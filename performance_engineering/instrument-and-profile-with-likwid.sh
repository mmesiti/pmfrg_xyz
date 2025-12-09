#!/bin/bash
julia --project=. -e 'using Pkg; Pkg.add("LIKWID"); Pkg.instantiate()'
# Apply instrumentation patch
git apply performance_engineering/likwid_markers.patch

# Run profiling (single-threaded, pinned to core 0)
cd performance_engineering
./run_likwid_profile_getXBubble.sh FLOPS_DP  

