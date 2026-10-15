#!/bin/bash
# All coarse runs (latitudes × grids × closures × seeds), single-threaded, 100 at a time.
# Create the job list first, then submit from the QGTestbed folder:
#   julia --project=. scripts/sweep.jl | grep -v truth | sed 's/^/JULIA_CPU_THREADS=1 /' > coarse_jobs.txt
#   sbatch scripts/slurm/coarse.sh
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
module load parallel 2>/dev/null   # GNU parallel, if not on the default path

# Compile QGTestbed once before starting the runs, so they don't all compile it at the same time
julia --project=. -e 'using QGTestbed' || { echo "QGTestbed does not load: run 'julia setup.jl' on the login node"; exit 1; }
export JULIA_NUM_THREADS=1

mkdir -p logs
parallel -j 100 --joblog coarse_joblog.txt --resume-failed --results logs/ < coarse_jobs.txt
