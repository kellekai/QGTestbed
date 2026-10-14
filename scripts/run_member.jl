# Run one sweep member from the command line, e.g.
#   julia --project=. scripts/run_member.jl 49 64 reynolds 20 5
#   julia --project=. scripts/run_member.jl 30 0 truth 20 5      # nx = 0 → auto truth resolution
# args: lat nx closure years spinup_years [seed] [gpu]
using QGTestbed
using GeophysicalFlows, FourierFlows

lat     = parse(Float64, ARGS[1])
nx      = parse(Int, ARGS[2])
closure = Symbol(ARGS[3])
years   = length(ARGS) ≥ 4 ? parse(Float64, ARGS[4]) : 20.0
spinup  = length(ARGS) ≥ 5 ? parse(Float64, ARGS[5]) : 5.0
seed    = length(ARGS) ≥ 6 ? parse(Int, ARGS[6]) : 1
use_gpu = length(ARGS) ≥ 7 && ARGS[7] == "gpu"
if use_gpu
    using CUDA
end
dev = use_gpu ? GPU() : CPU()

run_simulation(; lat, nx, closure, years, spinup_years = spinup, seed, dev,
               outdir = "output",
               # truth runs keep PV snapshots every 30 days for coarse-graining
               snapshot_every_days = closure === :truth ? 30.0 : Inf)
