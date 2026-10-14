# Write the full sweep as a list of commands (one per line) for a job array or
# GNU parallel:
#   julia --project=. scripts/sweep.jl > jobs.txt
#   parallel -j 4 < jobs.txt
#
# Sweep: latitudes × {truth, coarse grids × closures} × seeds.
lats     = 20:10:70
coarse   = (32, 64, 128)                 # 31, 15.6, 7.8 km for L = 1000 km
closures = (:none, :nemo, :geometric, :geometric_nosplit, :meke, :reynolds)
seeds    = 1:3                           # PG23 used 10 ensemble members
years, spinup = 20, 5
device = get(ENV, "QG_DEVICE", "cpu")   # set QG_DEVICE=gpu to run on CUDA

for lat in lats, seed in seeds
    println("julia --project=. scripts/run_member.jl $lat 0 truth $years $spinup $seed $device")
    for nx in coarse, cl in closures
        println("julia --project=. scripts/run_member.jl $lat $nx $cl $years $spinup $seed $device")
    end
end
