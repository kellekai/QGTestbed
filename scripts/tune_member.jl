# One run with given closure parameters (called from the job lists written by
# tune_plan.jl and scripts/slurm/tuned_sweep.sh):
#
#   julia --project=. scripts/tune_member.jl lat nx closure years spinup seed outdir [--no-tag] [NAME=value ...]
#
# The parameter values are stored in the output file (`metadata`) and, unless
# --no-tag is given, appended to its name, so runs with different settings never
# overwrite each other or the as-published sweep in output/. Without NAME=value
# pairs the closure runs with its published values.
using QGTestbed
include(joinpath(@__DIR__, "tune_params.jl"))

length(ARGS) >= 7 || error("usage: tune_member.jl lat nx closure years spinup seed outdir [--no-tag] [NAME=value ...]")
lat     = parse(Float64, ARGS[1])
nx      = parse(Int, ARGS[2])
closure = ARGS[3]
years   = parse(Float64, ARGS[4])
spinup  = parse(Float64, ARGS[5])
seed    = parse(Int, ARGS[6])
outdir  = ARGS[7]
rest    = ARGS[8:end]
notag   = "--no-tag" in rest
pairs   = [String(n) => parse(Float64, v) for (n, v) in (split(a, "="; limit = 2) for a in rest if a != "--no-tag")]

kw = isempty(pairs) ? NamedTuple() : closure_kwargs_from(closure, pairs)
run_simulation(; lat, nx, closure = Symbol(closure), years, spinup_years = spinup, seed, outdir,
               closure_kwargs = kw, tag = notag ? "" : tune_tag(pairs),
               metadata = Dict{String, Any}(n => v for (n, v) in pairs))
