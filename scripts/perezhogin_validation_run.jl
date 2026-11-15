# One run of the Perezhogin & Glazunov (2023) validation case (called by
# scripts/perezhogin_validation.sh):
#
#   julia --project=. scripts/perezhogin_validation_run.jl <member> <seed> <outdir> <years> <spinup> <C_R> \
#         [lat=49] [nx_truth=256] [nx_coarse=64] [snapshot_days=10] [aliased_fraction=0]
#
# aliased_fraction: 0 (default here) runs without GeophysicalFlows' 2/3 dealiasing, like pyqg,
# which relies on its exponential filter alone; 1/3 is the testbed default.
# member: truth    → no closure on nx_truth (the eddy-resolving reference)
#         none     → no closure on nx_coarse
#         reynolds → Reynolds-stress backscatter (Perezhogin & Glazunov 2023) with C_R on nx_coarse
# Full PV snapshots are stored every `snapshot_days` after the spin-up; the spectra and
# energy transfers are computed from them afterwards (perezhogin_validation_spectra.jl).
using QGTestbed

length(ARGS) >= 6 || error("usage: perezhogin_validation_run.jl member seed outdir years spinup C_R [lat] [nx_truth] [nx_coarse] [snapshot_days]")
member   = ARGS[1]
seed     = parse(Int, ARGS[2])
outdir   = ARGS[3]
years    = parse(Float64, ARGS[4])
spinup   = parse(Float64, ARGS[5])
C_R      = parse(Float64, ARGS[6])
lat      = length(ARGS) >= 7  ? parse(Float64, ARGS[7])  : 49.0
nx_truth = length(ARGS) >= 8  ? parse(Int, ARGS[8])      : 256
nx_c     = length(ARGS) >= 9  ? parse(Int, ARGS[9])      : 64
snapdays = length(ARGS) >= 10 ? parse(Float64, ARGS[10]) : 10.0
afrac    = length(ARGS) >= 11 ? parse(Float64, ARGS[11]) : 0.0   # 0 = no 2/3 dealiasing, as in pyqg

member in ("truth", "none", "reynolds") || error("member must be truth, none or reynolds, got '$member'")
nx = member == "truth" ? nx_truth : nx_c
closure_kwargs = member == "reynolds" ? (; C_R) : (;)

run_simulation(; lat, nx, closure = Symbol(member), years, spinup_years = spinup, seed, outdir,
               closure_kwargs, snapshot_every_days = snapdays, aliased_fraction = afrac,
               metadata = Dict{String, Any}("validation" => "Perezhogin & Glazunov 2023, Fig. 8",
                                            "C_R" => member == "reynolds" ? C_R : NaN))
