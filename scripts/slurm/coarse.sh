#!/bin/bash
# All coarse runs (latitudes × grids × closures × seeds), single-threaded, 100 at a time.
# Create the job list first, then submit from the QGTestbed folder:
#   julia --project=. scripts/sweep.jl | grep -v truth | sed 's/^/JULIA_CPU_THREADS=1 /' > coarse_jobs.txt
#   sbatch scripts/slurm/coarse.sh
# Rerun a subset, e.g. one closure after a code change (delete its output files first):
#   grep -E " (geometric|geometric_nosplit|meke) " coarse_jobs.txt > rerun_jobs.txt
#   JOBS=rerun_jobs.txt sbatch scripts/slurm/coarse.sh
# If the job hits the time limit, submit it again: --resume-failed skips successful runs and reruns failed ones.
#SBATCH --job-name=qg_coarse
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --time=48:00:00
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=coarse_%j.out

module load julia                  # the same Julia module used for `julia setup.jl`
# Numerics (README section 7d): no 2/3 dealiasing (validated against pyqg); exponential filter for
# all closures only if QG_FILTER_ALL=true (decide with scripts/slurm/stability.sh). Override when
# submitting, e.g.  QG_FILTER_ALL=true sbatch ...
export QG_ALIASED_FRACTION=${QG_ALIASED_FRACTION:-0}
export QG_FILTER_ALL=${QG_FILTER_ALL:-false}
echo "numerics: aliased_fraction = $QG_ALIASED_FRACTION, exponential filter for all closures = $QG_FILTER_ALL"
module load parallel 2>/dev/null   # GNU parallel, if not on the default path

# Compile QGTestbed once before starting the runs, so they don't all compile it at the same time
julia --project=. -e 'using QGTestbed' || { echo "QGTestbed does not load: run 'julia setup.jl' on the login node"; exit 1; }
export JULIA_NUM_THREADS=1

mkdir -p logs
# JOBS=<list> runs another job list (e.g. reruns of single closures) with its own job log;
# the default keeps the original names, so a resubmitted sweep continues where it stopped.
JOBS=${JOBS:-coarse_jobs.txt}
if [ "$JOBS" = coarse_jobs.txt ]; then JOBLOG=coarse_joblog.txt; else JOBLOG=${JOBS%.txt}_joblog.txt; fi
echo "job list $JOBS ($(wc -l < "$JOBS") runs), job log $JOBLOG"
parallel -j 100 --joblog "$JOBLOG" --resume-failed --results logs/ < "$JOBS"
