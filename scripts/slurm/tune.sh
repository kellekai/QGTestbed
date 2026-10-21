#!/bin/bash
# Tuning runs from a job list written by scripts/tune_plan.jl, single-threaded, 100 at a time.
#   julia scripts/tune_plan.jl --closure reynolds --param C_R=1,3,7,15,30 --cases 40:64,50:64 > tune_jobs.txt
#   sbatch scripts/slurm/tune.sh                 # or: JOBS=other_list.txt sbatch scripts/slurm/tune.sh
# Resubmit to finish interrupted lists: --resume-failed skips successful runs and reruns failed ones.
#SBATCH --job-name=qg_tune
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --time=48:00:00
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=tune_%j.out

module load julia                  # the same Julia module used for `julia setup.jl`
# Numerics (README section 7d): no 2/3 dealiasing (validated against pyqg); exponential filter for
# all closures only if QG_FILTER_ALL=true (decide with scripts/slurm/stability.sh). Override when
# submitting, e.g.  QG_FILTER_ALL=true sbatch ...
export QG_ALIASED_FRACTION=${QG_ALIASED_FRACTION:-0}
export QG_FILTER_ALL=${QG_FILTER_ALL:-false}
echo "numerics: aliased_fraction = $QG_ALIASED_FRACTION, exponential filter for all closures = $QG_FILTER_ALL"
module load parallel 2>/dev/null

# Compile QGTestbed once before starting the runs, so they don't all compile it at the same time
julia --project=. -e 'using QGTestbed' || { echo "QGTestbed does not load: run 'julia setup.jl' on the login node"; exit 1; }
export JULIA_NUM_THREADS=1

JOBS=${JOBS:-tune_jobs.txt}
mkdir -p logs/tune
parallel -j 100 --joblog ${JOBS%.txt}_joblog.txt --resume-failed --results logs/tune/ < "$JOBS"
