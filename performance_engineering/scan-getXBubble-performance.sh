#!/bin/bash
#SBATCH -n 1
#SBATCH --cpus-per-task=76
#SBATCH --partition cpuonly
#SBATCH --time 60
#SBATCH --cpu-freq=3000000

NTHREADS=76 # horeka

for commit in $(git log --format=%h ../src/PMFRG_xyz.jl)
do 
    git checkout $commit ../src/PMFRG_xyz.jl
    COMMENT="$(git log -1 $commit --format="%cn ; %s ; %ad")"
    julia --project=. --threads=$NTHREADS -O3 ./benchmark_and_test_getXBubble.jl --dbfile=getxbubble-scan.db --commit="$commit" --comment="$COMMENT"
done
