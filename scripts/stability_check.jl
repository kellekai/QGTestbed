# Stability check of the runs in one folder (written by scripts/slurm/stability.sh):
#
#   julia --project=plotting scripts/stability_check.jl -d output/stability_filter-false [-s 1]
#
# For every run:
#   status      ok | UNSTABLE (stopped on NaN, or KE not finite)
#   grid frac   fraction of the upper-layer KE above 0.65 k_max (k_max = π/Δx) in the last
#               vorticity snapshot: energy piling up at the grid scale
#   KE growth   KE1 at the end / median KE1 over the last 40 % of the run
# A run is flagged NOISY when its grid fraction exceeds both 0.02 and 3× that of the `none` run
# (exponential filter, no closure) at the same latitude and grid, and GROWING when its KE growth
# exceeds 3 and twice that of the `none` run (slow spin-up at high latitude, which the filtered
# reference shows too, is not counted). Writes <dir>/stability.csv and a per-closure summary.
ENV["GKSwstype"] = "100"
using JLD2, Printf, Statistics
include(joinpath(@__DIR__, "plotting_utils.jl"))     # YEAR, ke_spectrum_from_vorticity, closure_label

dir, spinup = "output/stability_filter-false", 1.0
let i = 1
    while i <= length(ARGS)
        a = ARGS[i]
        a in ("-d", "-s") && i < length(ARGS) || error("usage: stability_check.jl -d DIR [-s spinup_years]")
        a == "-d" && (global dir = ARGS[i+1])
        a == "-s" && (global spinup = parse(Float64, ARGS[i+1]))
        i += 2
    end
end
isdir(dir) || error("folder $dir does not exist")

const PATTERN = r"^(\d+)N_nx(\d+)_(.+)_seed(\d+)\.jld2$"
rows = []
for f in sort(readdir(dir))
    m = match(PATTERN, f)
    m === nothing && continue
    lat, nx, cl = parse(Int, m[1]), parse(Int, m[2]), m[3]
    d = load(joinpath(dir, f))
    ts = d["timeseries"]
    ke = ts["KE1"]; t = ts["t"] ./ YEAR
    stopped = get(d, "stopped_early", false) || any(!isfinite, ke)
    frac = NaN
    if !stopped && haskey(d, "vorticity_upper") && !isempty(d["vorticity_upper"])
        k, E = ke_spectrum_from_vorticity(d["vorticity_upper"][end], d["params"].L)
        kmax = π * nx / d["params"].L
        tot = sum(E)
        frac = tot > 0 ? sum(E[k .> 0.65 * kmax]) / tot : NaN
    end
    late = t .> max(spinup, 0.6 * t[end])
    growth = (!stopped && any(late)) ? ke[end] / median(ke[late]) : NaN
    push!(rows, (; lat, nx, closure = cl, stopped, frac, growth, stepper = get(d, "stepper", "?"),
                   afrac = get(d, "aliased_fraction", NaN)))
end
isempty(rows) && error("no run files in $dir")

ref = Dict((r.lat, r.nx) => r.frac for r in rows if r.closure == "none")
refg = Dict((r.lat, r.nx) => r.growth for r in rows if r.closure == "none")
function growing(r)
    (isfinite(r.growth) && r.growth > 3) || return false
    r.closure == "none" && return true
    g0 = get(refg, (r.lat, r.nx), NaN)
    return r.growth > 2 * (isfinite(g0) ? g0 : 1.0)
end
verdict(r) = r.stopped ? "UNSTABLE" :
             (isfinite(r.frac) && r.frac > max(0.02, 3 * get(ref, (r.lat, r.nx), 0.0))) ? "NOISY" :
             growing(r) ? "GROWING" : "ok"

open(joinpath(dir, "stability.csv"), "w") do io
    println(io, "lat,nx,closure,stepper,aliased_fraction,verdict,grid_fraction,ke_growth")
    for r in rows
        @printf(io, "%d,%d,%s,%s,%.4f,%s,%.4g,%.4g\n", r.lat, r.nx, r.closure, r.stepper, r.afrac, verdict(r), r.frac, r.growth)
    end
end

@printf("%-22s %5s %5s  %-12s %-9s %10s %10s\n", "closure", "lat", "grid", "stepper", "verdict", "grid frac", "KE growth")
for r in sort(rows; by = r -> (r.closure, r.lat, r.nx))
    @printf("%-22s %5d %5d  %-12s %-9s %10.4f %10.2f\n", closure_label(r.closure), r.lat, r.nx, r.stepper, verdict(r), r.frac, r.growth)
end
println("\nSummary (runs not ok / runs):")
for cl in sort(unique(r.closure for r in rows))
    rs = [r for r in rows if r.closure == cl]
    bad = [r for r in rs if verdict(r) != "ok"]
    @printf("  %-24s %d / %d %s\n", closure_label(cl), length(bad), length(rs),
            isempty(bad) ? "" : "(" * join(["$(r.lat)°N $(r.nx)² $(verdict(r))" for r in bad], ", ") * ")")
end
println("\nAll ok → keep QG_FILTER_ALL=false (closures as published). Otherwise rerun this job with")
println("QG_FILTER_ALL=true and, if that is clean, use QG_FILTER_ALL=true for the whole sweep.")
println("Saved $(joinpath(dir, "stability.csv"))")
