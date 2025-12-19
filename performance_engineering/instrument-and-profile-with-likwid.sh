#!/bin/bash
julia --project=. -e 'using Pkg; Pkg.add("LIKWID"); Pkg.instantiate()'
# Apply instrumentation patch
# git apply performance_engineering/likwid_markers.patch

# Run profiling (single-threaded, pinned to core 0)
cd performance_engineering


run_profiling(){
    PERFTYPE=$1
    PRECISION=$2
    echo "++++++++++++++++++++++++++++++++++++++++"
    echo "+ Checking FLOPS_${PERFTYPE}, executing FP${PRECISION}... +"
    echo "++++++++++++++++++++++++++++++++++++"
    ./run_likwid_profile_getXBubble.sh FLOPS_${PERFTYPE}  ${PRECISION}
}

run_profiling DP 64
run_profiling DP 32

run_profiling SP 64
run_profiling SP 32
