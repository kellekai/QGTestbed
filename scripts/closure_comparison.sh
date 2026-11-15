#!/bin/bash
# All closures (as published) at one latitude in the layout of Perezhogin & Glazunov (2023), Fig. 10:
# truth (256² at 40°N), and on 64² no closure, Reynolds averaging (NEMO default), Mak 2023,
# Jansen 2019 and Perezhogin 2023 (C_R = 7); then spectra, energy transfers and the figure.
# Submit from the QGTestbed root folder:
#   sbatch scripts/closure_comparison.sh
#   LAT=30 sbatch scripts/closure_comparison.sh                       # another latitude
#   CLOSURES="none nemo geometric geometric_nosplit meke reynolds" sbatch scripts/closure_comparison.sh
# Everything (runs, logs, spectra.jld2, closure_comparison.png, CSV curves) goes to $OUTDIR.
# Needs the model environment (julia setup.jl) and the plotting environment (julia setup.jl plotting).
# Rerunning skips runs whose output file already exists (e.g. after a time-out); delete them to start over.
# Reusing the sweep's truth runs (PV snapshots every 30 days, spin-up 5): link them into $OUTDIR and
# skip the truth runs; the spectra use only snapshots after $SPINUP years, as for the closures:
#   mkdir -p closure-comparison/40N_nx64 && ln -s $PWD/output/40N_nx256_truth_seed*.jld2 closure-comparison/40N_nx64/
#   TRUTH_SEEDS="" sbatch scripts/closure_comparison.sh
#SBATCH --job-name=qg_closure_cmp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --time=48:00:00
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=closure_comparison_%j.out

# ============================ settings ==========================================
LAT=${LAT:-40}
CLOSURES=${CLOSURES:-"none nemo nemo_eorca025 geometric meke reynolds"}   # nemo_eorca025 = our eORCA025 (no GM)
NX_COARSE=64
NX_TRUTH=0          # 0 = truth_resolution (4 points per Ld, at least 256²): 256² at 40°N
C_R=7               # Perezhogin & Glazunov (2023), published value
YEARS=20
SPINUP=8            # as in the proposal metrics (scale_awareness.jl -s 8)
SNAPSHOT_DAYS=30
SEEDS=${SEEDS:-"1 2 3 4 5 6 7 8 9 10"}
TRUTH_SEEDS=${TRUTH_SEEDS-$SEEDS}   # TRUTH_SEEDS="" runs no truth (use truth files linked into $OUTDIR)
TRUTH_THREADS=6     # 10 truth runs x 6 + 50 coarse runs x 1 = 110 cores
OUTDIR=${OUTDIR:-closure-comparison/${LAT}N_nx${NX_COARSE}}
# ================================================================================

module load julia                  # the same Julia module used for `julia setup.jl`
# Numerics of the sweep (README section 7d): no 2/3 dealiasing, exponential filter for every closure.
export QG_ALIASED_FRACTION=${QG_ALIASED_FRACTION:-0}
export QG_FILTER_ALL=${QG_FILTER_ALL:-true}
echo "numerics: aliased_fraction = $QG_ALIASED_FRACTION, exponential filter for all closures = $QG_FILTER_ALL"
set -o pipefail

julia --project=. -e 'using QGTestbed' || { echo "QGTestbed does not load: run 'julia setup.jl' on the login node"; exit 1; }
julia --project=plotting -e 'using JLD2, Plots' || { echo "plotting environment missing: run 'julia setup.jl plotting'"; exit 1; }

mkdir -p "$OUTDIR/logs"
{
    echo "Closure comparison, $(date)"
    echo "lat $LAT | closures $CLOSURES | coarse ${NX_COARSE}² | truth nx $NX_TRUTH (0 = automatic) | C_R $C_R | years $YEARS | spin-up $SPINUP | snapshots every $SNAPSHOT_DAYS d | seeds $SEEDS"
    echo "aliased_fraction $QG_ALIASED_FRACTION | filter all $QG_FILTER_ALL"
} > "$OUTDIR/settings.txt"
cat "$OUTDIR/settings.txt"

run() {   # member threads seed
    local member=$1 threads=$2 seed=$3
    local log="$OUTDIR/logs/${member}_seed${seed}.log"
    if ls "$OUTDIR"/${LAT}N_nx*_${member}_seed${seed}.jld2 >/dev/null 2>&1; then echo "exists, skipped: $member seed $seed"; return; fi
    JULIA_CPU_THREADS=$threads julia -t $threads --project=. scripts/closure_comparison_run.jl \
        $member $seed "$OUTDIR" $YEARS $SPINUP $LAT $NX_TRUTH $NX_COARSE $SNAPSHOT_DAYS $C_R > "$log" 2>&1
}

# 1) simulations, all in parallel
for seed in $SEEDS; do
    for cl in $CLOSURES; do run $cl 1 $seed & done
done
for seed in $TRUTH_SEEDS; do run truth $TRUTH_THREADS $seed & done
wait

missing=0
for seed in $SEEDS; do for m in $CLOSURES; do
    ls "$OUTDIR"/${LAT}N_nx*_${m}_seed${seed}.jld2 >/dev/null 2>&1 \
        || { echo "missing: $m seed $seed (see $OUTDIR/logs/${m}_seed${seed}.log)"; missing=1; }
done; done
ls "$OUTDIR"/${LAT}N_nx*_truth_seed*.jld2 >/dev/null 2>&1 || { echo "no truth run in $OUTDIR"; missing=1; }
[ $missing -eq 0 ] || echo "some runs failed; the analysis uses the runs that exist"

# 2) spectra and energy transfers (model environment)
JULIA_CPU_THREADS=8 julia -t 8 --project=. scripts/closure_comparison_spectra.jl "$OUTDIR" $NX_COARSE $SPINUP \
    2>&1 | tee "$OUTDIR/logs/spectra.log" || { echo "spectra failed"; exit 1; }

# 3) figures (plotting environment): all closures, and the three the proposal argues with
julia --project=plotting scripts/plot_closure_comparison.jl "$OUTDIR" 2>&1 | tee "$OUTDIR/logs/plot.log"
julia --project=plotting scripts/plot_closure_comparison.jl "$OUTDIR" none,nemo,reynolds 2>&1 | tee -a "$OUTDIR/logs/plot.log"
julia --project=plotting scripts/plot_closure_comparison.jl "$OUTDIR" none,nemo_eorca025,reynolds 2>&1 | tee -a "$OUTDIR/logs/plot.log"
echo "done: $OUTDIR/closure_comparison.png"
