#!/bin/bash
# Truth (eddy-resolving, no closure) runs for all latitudes, one seed, on one node.
# Submit from the QGTestbed folder:  sbatch scripts/slurm/truth.sh
#SBATCH --job-name=qg_truth
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --time=48:00:00
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=truth_%j.out

module load julia                  # the same Julia module used for `julia setup.jl`

# Compile QGTestbed once before starting the runs, so they don't all compile it at the same time
julia --project=. -e 'using QGTestbed' || { echo "QGTestbed does not load: run 'julia setup.jl' on the login node"; exit 1; }

SEED=${SEED:-1}
YEARS=${YEARS:-20}
SPINUP=${SPINUP:-5}
mkdir -p logs

# 512² grids (50–70°N): 16 FFT threads each. JULIA_CPU_THREADS sets the threads
# FourierFlows plans its FFTs for; -t gives Julia the same number.
for lat in 50 60 70; do
    JULIA_CPU_THREADS=16 julia -t 16 --project=. scripts/run_member.jl $lat 0 truth $YEARS $SPINUP $SEED \
        > logs/truth_${lat}N_seed${SEED}.log 2>&1 &
done
# 256² grids (20–40°N): 8 threads each.
for lat in 20 30 40; do
    JULIA_CPU_THREADS=8 julia -t 8 --project=. scripts/run_member.jl $lat 0 truth $YEARS $SPINUP $SEED \
        > logs/truth_${lat}N_seed${SEED}.log 2>&1 &
done
wait
