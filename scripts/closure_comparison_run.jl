# One run of the closure comparison in the layout of Perezhogin & Glazunov (2023), Fig. 10
# (called by scripts/closure_comparison.sh):
#
#   julia --project=. scripts/closure_comparison_run.jl <member> <seed> <outdir> <years> <spinup> \
#         [lat=40] [nx_truth=0] [nx_coarse=64] [snapshot_days=30] [C_R=7]
#
# member: truth                      → no closure on nx_truth (0 = truth_resolution, 256² at 40°N)
#         none | nemo | nemo_eorca025 | geometric | geometric_nosplit | geometric_eorca025 | meke | reynolds → that closure on nx_coarse,
#         settings as published (reynolds: C_R, default 7)
# Full PV snapshots (and, for geometric/meke, the eddy-energy field) are stored every
# `snapshot_days` after the spin-up; spectra and transfers are computed afterwards
# (closure_comparison_spectra.jl). Numerics follow QG_ALIASED_FRACTION / QG_FILTER_ALL.
using QGTestbed

length(ARGS) >= 5 || error("usage: closure_comparison_run.jl member seed outdir years spinup [lat] [nx_truth] [nx_coarse] [snapshot_days] [C_R]")
member   = ARGS[1]
seed     = parse(Int, ARGS[2])
outdir   = ARGS[3]
years    = parse(Float64, ARGS[4])
spinup   = parse(Float64, ARGS[5])
lat      = length(ARGS) >= 6  ? parse(Float64, ARGS[6])  : 40.0
nx_truth = length(ARGS) >= 7  ? parse(Int, ARGS[7])      : 0
nx_c     = length(ARGS) >= 8  ? parse(Int, ARGS[8])      : 64
snapdays = length(ARGS) >= 9  ? parse(Float64, ARGS[9])  : 30.0
C_R      = length(ARGS) >= 10 ? parse(Float64, ARGS[10]) : 7.0

members = ("truth", "none", "nemo", "nemo_eorca025", "geometric", "geometric_nosplit", "geometric_eorca025", "meke", "reynolds")
member in members || error("member must be one of $(join(members, ", ")), got '$member'")
nx = member == "truth" ? nx_truth : nx_c
closure_kwargs = member == "reynolds" ? (; C_R) : (;)

run_simulation(; lat, nx, closure = Symbol(member), years, spinup_years = spinup, seed, outdir,
               closure_kwargs, snapshot_every_days = snapdays,
               metadata = Dict{String, Any}("experiment" => "closure comparison, Perezhogin & Glazunov 2023 Fig. 10 layout"))
