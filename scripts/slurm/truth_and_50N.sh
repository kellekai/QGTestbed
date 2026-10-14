#!/bin/bash
# Truth runs for all latitudes plus the six 64² closure runs at 50°N (for a first
# vorticity comparison before the full coarse sweep), seed 1, on one node.
# Submit from the QGTestbed folder:  sbatch scripts/slurm/truth_and_50N.sh
# Afterwards remove the 50°N 64² seed-1 runs from coarse_jobs.txt (see README, section 5).
#SBATCH --job-name=qg_truth50
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --time=48:00:00
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=truth50_%j.out

module load julia                  # the same Julia module used for `julia setup.jl`
mkdir -p logs

# 512² truth runs: 16 threads each (48 cores)
for lat in 50 60 70; do
    JULIA_CPU_THREADS=16 julia -t 16 --project=. scripts/run_member.jl $lat 0 truth 20 5 1 > logs/truth_${lat}N.log 2>&1 &
done
# 256² truth runs: 8 threads each (24 cores)
for lat in 20 30 40; do
    JULIA_CPU_THREADS=8 julia -t 8 --project=. scripts/run_member.jl $lat 0 truth 20 5 1 > logs/truth_${lat}N.log 2>&1 &
done
# 64² closure runs at 50°N: 1 thread each (6 cores)
for cl in none nemo geometric geometric_nosplit meke reynolds; do
    JULIA_CPU_THREADS=1 julia -t 1 --project=. scripts/run_member.jl 50 64 $cl 20 5 1 > logs/${cl}_50N.log 2>&1 &
done

wait    # keep the job alive until every run above has finished
