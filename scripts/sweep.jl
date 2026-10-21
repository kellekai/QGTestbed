# Write the full sweep as a list of commands (one per line) for a job array or
# GNU parallel:
#   julia --project=. scripts/sweep.jl > jobs.txt
#   parallel -j 4 < jobs.txt
#
# Sweep: latitudes × {truth, coarse grids × closures} × seeds.
lats     = 20:10:70
coarse   = Tuple(parse.(Int, split(get(ENV, "QG_GRIDS", "64"), ",")))   # QG_GRIDS=32,64,128 for all (31, 15.6, 7.8 km)
closures = (:none, :nemo, :nemo_eorca025, :geometric, :geometric_nosplit, :geometric_eorca025, :meke, :reynolds)
seeds    = 1:3                           # PG23 used 10 ensemble members
years, spinup = 20, 5
device = get(ENV, "QG_DEVICE", "cpu")   # set QG_DEVICE=gpu to run on CUDA

for lat in lats, seed in seeds
    println("julia --project=. scripts/run_member.jl $lat 0 truth $years $spinup $seed $device")
    for nx in coarse, cl in closures
        println("julia --project=. scripts/run_member.jl $lat $nx $cl $years $spinup $seed $device")
    end
end
