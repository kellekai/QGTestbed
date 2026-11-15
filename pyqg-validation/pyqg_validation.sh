#!/bin/bash
# Perezhogin & Glazunov (2023) Fig. 10 in pyqg itself, with their own Reynolds_stress
# parameterization (parameterizations.py, same folder): 256² truth, 64² without and with the
# parameterization (C_R = 7), 10 members over years 5-20, then the same diagnostics as the
# QG testbed and the figure. Submit from this folder:
#   sbatch pyqg_validation.sh
# Needs Python with pyqg, numpy and matplotlib (pip install pyqg). Runs whose output file exists
# are skipped, so the job can be resubmitted after a time-out.
#SBATCH --job-name=pyqg_validation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --time=02:00:00
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=pyqg_validation_%j.out

# ============================ settings ==========================================
OUTDIR=results
SEEDS="1 2 3 4 5 6 7 8 9 10"
YEARS=20
SPINUP=5
SNAPSHOT_DAYS=30
CSIM=7
NX_TRUTH=256
NX_COARSE=64
TRUTH_THREADS=4      # pyfftw threads per truth run (ignored without pyfftw)
# Folder with the testbed's CSV curves to overlay (written by plot_perezhogin_validation.jl);
# leave empty for none, e.g. TESTBED=../perezhogin-validation/aliased_fraction_0
TESTBED=""
# ================================================================================

module load python 2>/dev/null     # or: source <your_venv>/bin/activate / conda activate <env>
set -o pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

python -c "import pyqg, numpy, matplotlib; print('pyqg', pyqg.__version__)" \
    || { echo "pyqg/numpy/matplotlib not importable: pip install pyqg matplotlib"; exit 1; }
[ -f parameterizations.py ] || { echo "parameterizations.py (Perezhogin, Zenodo 8277231) missing in $(pwd)"; exit 1; }

mkdir -p "$OUTDIR/logs"
run() {   # member nx threads seed
    local member=$1 nx=$2 threads=$3 seed=$4
    local file="$OUTDIR/49N_nx${nx}_${member}_seed${seed}.npz"
    if [ -f "$file" ]; then echo "exists, skipped: $file"; return; fi
    python run_pyqg.py --member $member --seed $seed --years $YEARS --spinup $SPINUP \
        --snapshot-days $SNAPSHOT_DAYS --nx-truth $NX_TRUTH --nx-coarse $NX_COARSE --csim $CSIM \
        --threads $threads --outdir "$OUTDIR" > "$OUTDIR/logs/${member}_seed${seed}.log" 2>&1
}

for seed in $SEEDS; do
    run truth    $NX_TRUTH  $TRUTH_THREADS $seed &
    run none     $NX_COARSE 1              $seed &
    run reynolds $NX_COARSE 1              $seed &
done
wait

ARGS=(--dir "$OUTDIR" --nx-coarse $NX_COARSE)
[ -n "$TESTBED" ] && ARGS+=(--testbed "$TESTBED")
python analyze_pyqg.py "${ARGS[@]}" 2>&1 | tee "$OUTDIR/logs/analyze.log"
echo "done: $OUTDIR/pyqg_validation.png"
