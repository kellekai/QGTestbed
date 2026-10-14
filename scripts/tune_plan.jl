# Write the job list for a tuning experiment: every combination of the given
# parameter values, for every case (latitude:grid) and seed.
#
#   julia scripts/tune_plan.jl --closure reynolds --param C_R=1,3,7,15,30 --cases 40:64,50:64 > tune_jobs.txt
#   julia scripts/tune_plan.jl --closure geometric --param alpha=0.02,0.04,0.06,0.1 --param lambda_days=20,80,320 \
#         --cases 40:64 --years 10 --spinup 3 > tune_jobs.txt
#
# Options
#   --closure NAME        nemo | geometric | geometric_nosplit | meke | reynolds
#   --param NAME=v1,v2,…  values to try (repeat --param for more parameters → full factorial);
#                         names are listed in scripts/tune_params.jl
#   --cases L:N,…         latitude:grid pairs to run, e.g. 40:64,60:32 (default 40:64)
#   --years Y             length of each run (default 10)
#   --spinup S            spin-up excluded from the scores (default 3); energies level off after ~2 yr
#   --seeds 1,2           seeds (default 1)
#   --outdir DIR          where runs are written (default output/tune/<closure>)
#
# Run the job list like the coarse sweep (scripts/slurm/tune.sh, GNU parallel),
# then score it with scripts/tune_evaluate.jl. Needs no packages.
include(joinpath(@__DIR__, "tune_params.jl"))

closure, cases, years, spinup, seeds, outdir = "", ["40:64"], 10.0, 3.0, [1], ""
params = Pair{String, Vector{Float64}}[]
let i = 1
    while i <= length(ARGS)
        a = ARGS[i]
        i < length(ARGS) || error("option $a needs a value")
        v = ARGS[i+1]
        if a == "--closure"
            global closure = v
        elseif a == "--param"
            name, vals = split(v, "="; limit = 2)
            push!(params, String(name) => parse.(Float64, split(vals, ",")))
        elseif a == "--cases"
            global cases = String.(split(v, ","))
        elseif a == "--years"
            global years = parse(Float64, v)
        elseif a == "--spinup"
            global spinup = parse(Float64, v)
        elseif a == "--seeds"
            global seeds = parse.(Int, split(v, ","))
        elseif a == "--outdir"
            global outdir = v
        else
            error("unknown option $a")
        end
        i += 2
    end
end
isempty(closure) && error("--closure is required")
isempty(params) && error("give at least one --param NAME=v1,v2,...")
outdir = isempty(outdir) ? joinpath("output", "tune", closure) : outdir

# validate names (and the conversions) once, before writing anything
closure_kwargs_from(closure, [n => first(v) for (n, v) in params])
for c in cases
    occursin(r"^\d+(\.\d+)?:\d+$", c) || error("case '$c' must be latitude:grid, e.g. 40:64")
end

combos = vec(collect(Iterators.product((v for (_, v) in params)...)))
pnames = first.(params)
njobs = 0
for c in cases, combo in combos, seed in seeds
    lat, nx = split(c, ":")
    settings = join(["$(n)=$(v)" for (n, v) in zip(pnames, combo)], " ")
    println("JULIA_CPU_THREADS=1 julia -t 1 --project=. scripts/tune_member.jl $lat $nx $closure $years $spinup $seed $outdir $settings")
    global njobs += 1
end
println(stderr, "wrote $njobs jobs: $(length(combos)) parameter sets × $(length(cases)) cases × $(length(seeds)) seeds → $outdir")
