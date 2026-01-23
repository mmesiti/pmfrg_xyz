#!/usr/bin/env sh

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4 # FIX THIS
#SBATCH --mem-per-cpu=500M

#SBATCH --output=./logging/slurm_%j.out
#SBATCH --error=./logging/slurm_%j.err
#SBATCH --time=0-02:00:00
#SBATCH --partition=main

set -eu
#scontrol show job "$SLURM_JOBID"

#module purge
#module load julia

JULIA_SCRIPT="$1"
julia --project=. -t "${SLURM_CPUS_PER_TASK:-4}" "$JULIA_SCRIPT"

