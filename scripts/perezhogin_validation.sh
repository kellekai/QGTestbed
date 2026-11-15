#!/bin/bash
# Validation against Perezhogin & Glazunov (2023), Fig. 10 (Fig. 8 in the preprint): the pyqg "eddy" configuration
# (49°N member of the testbed), truth on 256², 64² without closure and 64² with the
# Reynolds-stress backscatter (C_R = 7), then spectra, energy transfers and the figure.
# Submit from the QGTestbed root folder:
#   sbatch scripts/perezhogin_validation.sh
# Everything (runs, logs, spectra.jld2, perezhogin_validation.png) goes to $OUTDIR, a subfolder of
# perezhogin-validation/ named after the dealiasing setting, so runs with other settings stay apart.
# Needs the model environment (julia setup.jl) and the plotting environment (julia setup.jl plotting).
# Rerunning skips runs whose output file already exists; delete them to start over.
#SBATCH --job-name=qg_perezhogin_val
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --time=24:00:00
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=perezhogin_validation_%j.out

# ============================ settings ==========================================
ALIASED_FRACTION=0  # 0 = no 2/3 dealiasing, exponential filter only, as in pyqg; 0.3333333 = testbed default
OUTDIR=perezhogin-validation/aliased_fraction_${ALIASED_FRACTION}
LAT=49            # 49°N = the pyqg eddy configuration (Ld = 15 km, beta = 1.5e-11, dU = 0.025 m/s)
NX_TRUTH=256      # as in Ross et al. (2023) / Perezhogin & Glazunov (2023)
NX_COARSE=64
C_R=7             # published value for this setup
YEARS=20          # Perezhogin & Glazunov (2023): averages over years 5-20
SPINUP=5
SNAPSHOT_DAYS=30  # PV snapshots for the spectra, after the spin-up (181 per run)
SEEDS="1 2 3 4 5 6 7 8 9 10"   # ensemble of 10 members, as in the paper
TRUTH_THREADS=8   # 10 truth runs x 8 + 20 coarse runs x 1 = 100 cores
# ================================================================================

module load julia                  # the same Julia module used for `julia setup.jl`
set -o pipefail                    # a failing julia step is not hidden by `tee`

julia --project=. -e 'using QGTestbed' || { echo "QGTestbed does not load: run 'julia setup.jl' on the login node"; exit 1; }
julia --project=plotting -e 'using JLD2, Plots' || { echo "plotting environment missing: run 'julia setup.jl plotting'"; exit 1; }

mkdir -p "$OUTDIR/logs"
{
    echo "Perezhogin & Glazunov (2023) validation, $(date)"
    echo "lat $LAT | truth ${NX_TRUTH}² | coarse ${NX_COARSE}² | C_R $C_R | years $YEARS | spin-up $SPINUP | snapshots every $SNAPSHOT_DAYS d | seeds $SEEDS | aliased_fraction $ALIASED_FRACTION"
} > "$OUTDIR/settings.txt"
cat "$OUTDIR/settings.txt"

run() {   # member nx threads seed
    local member=$1 nx=$2 threads=$3 seed=$4
    local file="$OUTDIR/${LAT}N_nx${nx}_${member}_seed${seed}.jld2"
    if [ -f "$file" ]; then echo "exists, skipped: $file"; return; fi
    JULIA_CPU_THREADS=$threads julia -t $threads --project=. scripts/perezhogin_validation_run.jl \
        $member $seed "$OUTDIR" $YEARS $SPINUP $C_R $LAT $NX_TRUTH $NX_COARSE $SNAPSHOT_DAYS $ALIASED_FRACTION \
        > "$OUTDIR/logs/${member}_nx${nx}_seed${seed}.log" 2>&1
}

# 1) simulations, all in parallel
for seed in $SEEDS; do
    run truth    $NX_TRUTH  $TRUTH_THREADS $seed &
    run none     $NX_COARSE 1              $seed &
    run reynolds $NX_COARSE 1              $seed &
done
wait

missing=0
for seed in $SEEDS; do for m in "truth $NX_TRUTH" "none $NX_COARSE" "reynolds $NX_COARSE"; do
    set -- $m
    f="$OUTDIR/${LAT}N_nx${2}_${1}_seed${seed}.jld2"
    [ -f "$f" ] || { echo "missing $f (see $OUTDIR/logs/${1}_nx${2}_seed${seed}.log)"; missing=1; }
done; done
[ $missing -eq 0 ] || echo "some runs failed; the analysis uses the runs that exist"

# 2) spectra and energy transfers (model environment: needs the closure and the coarse-graining)
JULIA_CPU_THREADS=8 julia -t 8 --project=. scripts/perezhogin_validation_spectra.jl "$OUTDIR" $NX_COARSE \
    2>&1 | tee "$OUTDIR/logs/spectra.log" || { echo "spectra failed"; exit 1; }

# 3) figure (plotting environment)
julia --project=plotting scripts/plot_perezhogin_validation.jl "$OUTDIR" 2>&1 | tee "$OUTDIR/logs/plot.log"
echo "done: $OUTDIR/perezhogin_validation.png"
