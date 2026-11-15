#!/bin/bash
# Repeat the coarse sweep with tuned closure parameters.
# Edit the settings block below (best values from tune_evaluate.jl), then submit from the
# QGTestbed folder:   sbatch scripts/slurm/tuned_sweep.sh
#
# Runs go to $OUTDIR (default output/tuned/) with the same file names as the as-published
# sweep in output/, which stays untouched. The parameter values are stored in every file
# (`metadata`) and in $OUTDIR/tuned_settings.txt. Compare with the truth runs in output/:
#   julia --project=plotting scripts/scale_awareness.jl -d output/tuned -t output
# Resubmit with the same settings to finish an interrupted sweep (--resume-failed).
#SBATCH --job-name=qg_tuned_sweep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --time=48:00:00
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=tuned_sweep_%j.out

# ============================ settings ==========================================
# Best parameter values per closure, space-separated NAME=value (names as in tune_plan.jl:
#   nemo: Uv Ue Le_km | geometric(_nosplit): alpha lambda_days etaE E0 L_split_km Uv
#   meke: c c_smag tau_vis_days c_d U_bg | reynolds: C_R FGR).
# Parameters not listed keep their published values. An empty string runs the closure
# with its published values; remove a closure from CLOSURES to skip it.
NEMO_PARAMS="Uv=0.1 Ue=0.02"
GEOMETRIC_PARAMS="alpha=0.06 lambda_days=80"
GEOMETRIC_NOSPLIT_PARAMS="alpha=0.06 lambda_days=80"
MEKE_PARAMS="c=0.15"
REYNOLDS_PARAMS="C_R=3"

CLOSURES="none nemo geometric geometric_nosplit meke reynolds"   # `none` = reference without closure
LATS="20 30 40 50 60 70"
GRIDS="64"                       # the grid used for the proposal; "32 64 128" for all
SEEDS="1 2 3"
YEARS=20
SPINUP=5
OUTDIR=output/tuned
NPAR=100                                                          # runs at a time
# ================================================================================

module load julia                  # the same Julia module used for `julia setup.jl`
# Numerics (README section 7d): no 2/3 dealiasing (validated against pyqg); exponential filter for
# all closures only if QG_FILTER_ALL=true (decide with scripts/slurm/stability.sh). Override when
# submitting, e.g.  QG_FILTER_ALL=true sbatch ...
export QG_ALIASED_FRACTION=${QG_ALIASED_FRACTION:-0}
export QG_FILTER_ALL=${QG_FILTER_ALL:-false}
echo "numerics: aliased_fraction = $QG_ALIASED_FRACTION, exponential filter for all closures = $QG_FILTER_ALL"
module load parallel 2>/dev/null

params_for() {
    case "$1" in
        nemo)              echo "$NEMO_PARAMS" ;;
        geometric)         echo "$GEOMETRIC_PARAMS" ;;
        geometric_nosplit) echo "$GEOMETRIC_NOSPLIT_PARAMS" ;;
        meke)              echo "$MEKE_PARAMS" ;;
        reynolds)          echo "$REYNOLDS_PARAMS" ;;
        *)                 echo "" ;;
    esac
}

# Compile QGTestbed once before starting the runs, so they don't all compile it at the same time
julia --project=. -e 'using QGTestbed' || { echo "QGTestbed does not load: run 'julia setup.jl' on the login node"; exit 1; }
export JULIA_NUM_THREADS=1

# Check the parameter names before starting hundreds of runs (tune_plan.jl validates them)
for cl in $CLOSURES; do
    P=$(params_for $cl)
    [ -z "$P" ] && continue
    ARGS=(); for p in $P; do ARGS+=(--param "$p"); done
    julia scripts/tune_plan.jl --closure $cl "${ARGS[@]}" > /dev/null 2>&1 \
        || { echo "invalid parameters for $cl: '$P' (see scripts/tune_params.jl)"; exit 1; }
done

# Record the settings and write the job list
mkdir -p "$OUTDIR" logs/tuned
{
    echo "tuned sweep, $(date)"
    for cl in $CLOSURES; do
        P=$(params_for $cl); echo "$cl: ${P:-published values}"
    done
    echo "lats: $LATS | grids: $GRIDS | seeds: $SEEDS | years: $YEARS | spin-up: $SPINUP"
} > "$OUTDIR/tuned_settings.txt"
cat "$OUTDIR/tuned_settings.txt"

JOBS=tuned_sweep_jobs.txt
> "$JOBS"
for lat in $LATS; do for nx in $GRIDS; do for seed in $SEEDS; do for cl in $CLOSURES; do
    echo "JULIA_CPU_THREADS=1 julia -t 1 --project=. scripts/tune_member.jl $lat $nx $cl $YEARS $SPINUP $seed $OUTDIR --no-tag $(params_for $cl)" >> "$JOBS"
done; done; done; done
echo "$(wc -l < "$JOBS") runs in $JOBS"

parallel -j $NPAR --joblog tuned_sweep_joblog.txt --resume-failed --results logs/tuned/ < "$JOBS"
