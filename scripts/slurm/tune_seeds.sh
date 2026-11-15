#!/bin/bash
# Tuning runs with several seeds: writes its own job list with scripts/tune_plan.jl and runs it.
# All settings come from environment variables, so no job list has to be prepared beforehand.
#
#   CLOSURE=reynolds PARAMS="C_R=0,1,3,7,15,30" CASES=40:64,60:64 SEEDS=2,3 \
#       sbatch scripts/slurm/tune_seeds.sh
#
#   CLOSURE=geometric PARAMS="alpha=0.02,0.04,0.06,0.1 lambda_days=20,80,320" \
#       CASES=40:64,60:64 SEEDS=1,2,3 sbatch scripts/slurm/tune_seeds.sh
#
#   CLOSURE   closure to tune (nemo | geometric | geometric_nosplit | meke | reynolds)   [required]
#   PARAMS    space-separated NAME=v1,v2,… (one per parameter, all combinations run)    [required]
#   CASES     latitude:grid pairs                                         (default 40:64,60:64)
#   SEEDS     seeds to run                                                (default 1,2,3)
#   YEARS     run length in years                                         (default 10)
#   SPINUP    spin-up excluded from the scores                            (default 3)
#
# Runs go to output/tune/<closure>/ like the single-seed tuning, so seeds already run there
# (e.g. seed 1 from scripts/slurm/tune.sh) are reused: set SEEDS=2,3 to add only new seeds.
# tune_evaluate.jl averages the cost over all seeds it finds.
# Resubmit with the same settings to finish an interrupted list (--resume-failed).
#SBATCH --job-name=qg_tune_seeds
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --time=48:00:00
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=tune_seeds_%j.out

module load julia                  # the same Julia module used for `julia setup.jl`
# Numerics (README section 7d): no 2/3 dealiasing (validated against pyqg); exponential filter for
# all closures only if QG_FILTER_ALL=true (decide with scripts/slurm/stability.sh). Override when
# submitting, e.g.  QG_FILTER_ALL=true sbatch ...
export QG_ALIASED_FRACTION=${QG_ALIASED_FRACTION:-0}
export QG_FILTER_ALL=${QG_FILTER_ALL:-false}
echo "numerics: aliased_fraction = $QG_ALIASED_FRACTION, exponential filter for all closures = $QG_FILTER_ALL"
module load parallel 2>/dev/null

: "${CLOSURE:?set CLOSURE, e.g. CLOSURE=reynolds}"
: "${PARAMS:?set PARAMS, e.g. PARAMS=\"C_R=0,1,3,7,15,30\"}"
CASES=${CASES:-40:64,60:64}
SEEDS=${SEEDS:-1,2,3}
YEARS=${YEARS:-10}
SPINUP=${SPINUP:-3}

# Compile QGTestbed once before starting the runs, so they don't all compile it at the same time
julia --project=. -e 'using QGTestbed' || { echo "QGTestbed does not load: run 'julia setup.jl' on the login node"; exit 1; }
export JULIA_NUM_THREADS=1

# Job list for these settings (one file per closure and seed set)
JOBS=tune_${CLOSURE}_seeds_${SEEDS//,/-}.txt
PARAM_ARGS=()
for p in $PARAMS; do PARAM_ARGS+=(--param "$p"); done
julia scripts/tune_plan.jl --closure "$CLOSURE" "${PARAM_ARGS[@]}" --cases "$CASES" \
      --seeds "$SEEDS" --years "$YEARS" --spinup "$SPINUP" > "$JOBS" || { echo "tune_plan.jl failed"; exit 1; }
echo "$(wc -l < "$JOBS") runs in $JOBS"

mkdir -p logs/tune
parallel -j 100 --joblog ${JOBS%.txt}_joblog.txt --resume-failed --results logs/tune/ < "$JOBS"
