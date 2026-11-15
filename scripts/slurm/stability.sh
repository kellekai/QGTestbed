#!/bin/bash
# Stability test of the closures without 2/3 dealiasing (README section 7d).
# nemo, geometric, geometric_nosplit and meke run with plain RK4 and their own dissipation; with
# dealiasing switched off nothing else removes grid-scale aliasing. This job runs every closure
# (published settings) for a few years on 64² at all six latitudes (LATS/GRIDS to change), then
# scripts/stability_check.jl reports runs that blew up (NaN) or pile up energy at the grid scale.
# Submit from the QGTestbed folder:
#   sbatch scripts/slurm/stability.sh                       # closures as published (plain RK4)
#   QG_FILTER_ALL=true sbatch scripts/slurm/stability.sh    # every closure with the exponential filter
# Results go to output/stability_filter-<QG_FILTER_ALL>/, a table to stability.csv in that folder.
#SBATCH --job-name=qg_stability
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --time=06:00:00
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=stability_%j.out

# ============================ settings ==========================================
LATS=${LATS:-"20 30 40 50 60 70"}   # override when submitting, e.g. GRIDS="32 64 128" sbatch ...
GRIDS=${GRIDS:-"64"}                # the grid used for the proposal
CLOSURES=${CLOSURES:-"none nemo geometric geometric_nosplit meke reynolds"}   # none = filtered reference
YEARS=5
SPINUP=1
SEED=1
NPAR=100
# ================================================================================

module load julia                  # the same Julia module used for `julia setup.jl`
export QG_ALIASED_FRACTION=${QG_ALIASED_FRACTION:-0}
export QG_FILTER_ALL=${QG_FILTER_ALL:-false}
echo "numerics: aliased_fraction = $QG_ALIASED_FRACTION, exponential filter for all closures = $QG_FILTER_ALL"
module load parallel 2>/dev/null
set -o pipefail

julia --project=. -e 'using QGTestbed' || { echo "QGTestbed does not load: run 'julia setup.jl' on the login node"; exit 1; }
export JULIA_NUM_THREADS=1

OUTDIR=${OUTDIR:-output/stability_filter-${QG_FILTER_ALL}}
mkdir -p "$OUTDIR" logs/stability
JOBS=$OUTDIR/stability_jobs.txt
> "$JOBS"
for lat in $LATS; do for nx in $GRIDS; do for cl in $CLOSURES; do
    echo "JULIA_CPU_THREADS=1 julia -t 1 --project=. scripts/tune_member.jl $lat $nx $cl $YEARS $SPINUP $SEED $OUTDIR --no-tag" >> "$JOBS"
done; done; done
echo "$(wc -l < "$JOBS") runs"

parallel -j $NPAR --joblog "$OUTDIR/stability_joblog.txt" --resume-failed --results logs/stability/ < "$JOBS"

julia --project=plotting scripts/stability_check.jl -d "$OUTDIR" -s $SPINUP 2>&1 | tee "$OUTDIR/stability_summary.txt"
