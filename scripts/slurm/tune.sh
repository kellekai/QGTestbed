#!/bin/bash
# Tuning runs from a job list written by scripts/tune_plan.jl, single-threaded, 100 at a time.
#   julia scripts/tune_plan.jl --closure reynolds --param C_R=1,3,7,15,30 --cases 40:64,50:64 > tune_jobs.txt
#   sbatch scripts/slurm/tune.sh                 # or: JOBS=other_list.txt sbatch scripts/slurm/tune.sh
# Resubmit to finish interrupted lists: --resume skips completed runs.
#SBATCH --job-name=qg_tune
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --time=48:00:00
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=tune_%j.out

module load julia                  # the same Julia module used for `julia setup.jl`
module load parallel 2>/dev/null
export JULIA_NUM_THREADS=1

JOBS=${JOBS:-tune_jobs.txt}
mkdir -p logs/tune
parallel -j 100 --joblog ${JOBS%.txt}_joblog.txt --resume --results logs/tune/ < "$JOBS"
