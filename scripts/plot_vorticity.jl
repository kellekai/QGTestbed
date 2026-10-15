# Upper-layer relative vorticity ζ₁/f₀ for every closure at one latitude and
# grid, side by side with the coarse-grained truth.
#
#   julia --project=plotting scripts/plot_vorticity.jl 50 64              # times chosen automatically
#   julia --project=plotting scripts/plot_vorticity.jl 50 64 0.5 1 2 5    # times in years
#   julia --project=plotting scripts/plot_vorticity.jl 50 64 -d output/tuned -t output -o figures
#
# Options (anywhere after the script name):
#   -d DIR   folder with the coarse runs (default output)
#   -t DIR   folder with the truth runs (default: same as -d); the no-closure runs are
#            also taken from here if the -d folder has none
#   -o DIR   folder for the figure (default: same as -d; created if missing)
#
# Uses the separate plotting environment (create it once with `julia setup.jl plotting`);
# it does not load QGTestbed, only the .jld2 output files. Plots.jl with the GR
# backend, which runs headless on cluster nodes.
#
# Columns: truth (coarse-grained to the same grid), then the closures, titled with
# the names in LABELS (scripts/plotting_utils.jl); missing runs are skipped.
# Rows: times. Automatic times: t* = first time the upper-layer KE of the
# reference run (truth if present, otherwise `none`) reaches half its maximum,
# then t* + 1 yr, t* + 3 yr and the last snapshot.
# One colour scale per row, symmetric, set by the 99.5th percentile of |ζ/f₀|
# of the reference panel in that row, so damped closures appear faded.
# Output: <-o DIR>/vorticity_<lat>N_nx<nx>.png
ENV["GKSwstype"] = "100"            # GR without a display
using JLD2, Plots, FFTW, Statistics, Printf
gr()

include(joinpath(@__DIR__, "plotting_utils.jl"))   # YEAR, coarsegrain2d, closure_label
const CLOSURES = ["none", "nemo", "geometric", "geometric_nosplit", "meke", "reynolds"]

pos = String[]; outdir = "output"; truthdir = ""; figdir = ""
let i = 1
    while i <= length(ARGS)
        a = ARGS[i]
        if a in ("-d", "-t", "-o")
            i < length(ARGS) || error("option $a needs a folder")
            a == "-d" && (global outdir = ARGS[i+1])
            a == "-t" && (global truthdir = ARGS[i+1])
            a == "-o" && (global figdir = ARGS[i+1])
            i += 2
        else
            push!(pos, a); i += 1
        end
    end
end
length(pos) >= 2 || error("usage: plot_vorticity.jl <lat> <nx> [times in years...] [-d DIR] [-t DIR] [-o DIR]")
lat = parse(Float64, pos[1])
nx  = parse(Int, pos[2])
times_requested = length(pos) > 2 ? parse.(Float64, pos[3:end]) : Float64[]
truthdir = isempty(truthdir) ? outdir : truthdir
figdir = isempty(figdir) ? outdir : figdir
tag = @sprintf("%02.0fN", lat)
@info "Arguments" latitude = lat grid = nx times = (isempty(times_requested) ? "automatic" : times_requested)
let available = sort(unique(parse(Int, m[1]) for f in readdir(outdir)
                            for m in (match(Regex("^$(tag)_nx(\\d+)_(?!truth).+_seed1\\.jld2\$"), f),) if m !== nothing))
    nx in available || error("no coarse runs for $tag on a $(nx)² grid in $outdir; available grids for $tag: $available")
end

function find_truth(outdir, tag)
    files = filter(f -> startswith(f, tag * "_nx") && endswith(f, "_truth_seed1.jld2"), readdir(outdir))
    return isempty(files) ? nothing : joinpath(outdir, first(files))
end


# --- load ----------------------------------------------------------------------
columns = Tuple{String, Dict{String, Any}}[]
truthfile = find_truth(truthdir, tag)
truthfile === nothing || push!(columns, ("truth (coarse-grained)", load(truthfile)))
for cl in CLOSURES
    file = joinpath(outdir, @sprintf("%s_nx%d_%s_seed1.jld2", tag, nx, cl))
    if !isfile(file) && cl == "none" && truthdir != outdir     # reference from the truth folder
        file = joinpath(truthdir, @sprintf("%s_nx%d_%s_seed1.jld2", tag, nx, cl))
    end
    isfile(file) ? push!(columns, (closure_label(cl), load(file))) : @warn "missing $file, column skipped"
end
isempty(columns) && error("no output files for $tag at nx = $nx in $outdir")
for (name, d) in columns
    haskey(d, "vorticity_upper") ||
        error("run '$name' has no vorticity snapshots: it was made with an older run_simulation, rerun it")
end

ref = columns[1][2]
pp  = ref["params"]
f0, L = pp.f₀, pp.L

# --- times -----------------------------------------------------------------------
function auto_times(d)
    t, KE = d["timeseries"]["t"] ./ YEAR, d["timeseries"]["KE1"]
    i = findfirst(>=(0.5 * maximum(KE)), KE)
    tstar = t[i]
    tend = d["vorticity_times"][end] / YEAR
    return unique(filter(<=(tend), [tstar, tstar + 1, tstar + 3, tend]))
end
times = isempty(times_requested) ? auto_times(ref) : times_requested
@info "Plotting times (years)" times

nearest(d, t) = argmin(abs.(d["vorticity_times"] ./ YEAR .- t))

# --- figure ----------------------------------------------------------------------
ncol, nrow = length(columns), length(times)
xkm = collect(range(-L / 2, step = L / nx, length = nx)) ./ 1e3
panels = Plots.Plot[]

for (r, t) in enumerate(times)
    fields = [coarsegrain2d(d["vorticity_upper"][nearest(d, t)], nx, L) ./ f0 for (_, d) in columns]
    a = quantile(vec(abs.(fields[1])), 0.995)
    a = a > 0 ? a : 1e-12
    for (c, (name, d)) in enumerate(columns)
        tshow = d["vorticity_times"][nearest(d, t)] / YEAR
        p = heatmap(xkm, xkm, permutedims(fields[c]);          # Plots wants z[y, x]
                    c = :balance, clims = (-a, a), aspect_ratio = :equal,
                    xlims = extrema(xkm), ylims = extrema(xkm), colorbar = false,
                    title = r == 1 ? name : "", titlefontsize = 9, tickfontsize = 6,
                    xticks = r == nrow ? :auto : false, yticks = c == 1 ? :auto : false,
                    xlabel = r == nrow ? "x [km]" : "", ylabel = c == 1 ? "y [km]" : "",
                    guidefontsize = 7, framestyle = :box)
        annotate!(p, xkm[1] + 0.04 * (xkm[end] - xkm[1]), xkm[end] - 0.08 * (xkm[end] - xkm[1]),
                  text(@sprintf("%.1f yr", tshow), 7, :left))
        push!(panels, p)
    end
    # Colour scale of this row as its own narrow panel (a vertical gradient strip
    # with value ticks), so every vorticity panel has the same size.
    levels = collect(range(-a, a, length = 201))
    push!(panels, heatmap([0.0, 1.0], levels, repeat(levels, 1, 2);
                          c = :balance, clims = (-a, a), colorbar = false,
                          xticks = false, ymirror = true, ylims = (-a, a),
                          tickfontsize = 6, framestyle = :box, right_margin = 8Plots.mm,
                          title = r == 1 ? "ζ₁/f₀" : "", titlefontsize = 9))
end

wpanel = 0.95 / ncol
fig = plot(panels...; layout = grid(nrow, ncol + 1, widths = [fill(wpanel, ncol); 1 - ncol * wpanel]),
           size = (270 * ncol + 90, 250 * nrow + 40),
           plot_title = @sprintf("Upper-layer relative vorticity, %.0f°N, %d² (Δx/Ld = %.2f)",
                                 lat, nx, (L / nx) / pp.Ld),
           plot_titlefontsize = 11)

mkpath(figdir)
file = joinpath(figdir, @sprintf("vorticity_%s_nx%d.png", tag, nx))
savefig(fig, file)
@info "Saved $file"
