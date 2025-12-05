#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=76
#SBATCH --partition cpuonly
#SBATCH --time 25
#SBATCH --cpu-freq=3000000

set -eu

julia --threads=$SLURM_CPUS_PER_TASK -O3 ./profile.jl square_lattice_large
