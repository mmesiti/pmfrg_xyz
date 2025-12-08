#!/bin/bash
# Run LIKWID profiling for PMFRG_xyz
# Usage: ./run_likwid_profile.sh [GROUP]
# Example: ./run_likwid_profile.sh FLOPS_DP

set -e

main() {
    local group="${1:-FLOPS_DP}"
    local script_dir
    script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    local project_root
    project_root="$(dirname "$script_dir")"

    cd "$script_dir"

    echo "=== LIKWID Profiling ==="
    echo "Group: $group"
    echo "Project: $project_root"
    echo ""

    likwid-perfctr -C 0 -g "$group" -m \
        julia --project="$project_root" -O3 -t 1 likwid_profile_getXBubble.jl "$group"
}

main "$@"
