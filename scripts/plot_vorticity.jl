# Upper-layer relative vorticity ζ₁/f₀ for every closure at one latitude and
# grid, side by side with the coarse-grained truth.
#
#   julia --project=plotting scripts/plot_vorticity.jl 50 64              # times chosen automatically
#   julia --project=plotting scripts/plot_vorticity.jl 50 64 0.5 1 2 5    # times in years
#
# Uses the separate plotting environment (create it once with `julia setup.jl plotting`);
# it does not load QGTestbed, only the .jld2 output files. Plots.jl with the GR
# backend, which runs headless on cluster nodes.
#
# Columns: truth (coarse-grained to the same grid), none, nemo, geometric,
# geometric_nosplit, meke, reynolds (missing runs are skipped).
# Rows: times. Automatic times: t* = first time the upper-layer KE of the
# reference run (truth if present, otherwise `none`) reaches half its maximum,
# then t* + 1 yr, t* + 3 yr and the last snapshot.
# One colour scale per row, symmetric, set by the 99.5th percentile of |ζ/f₀|
# of the reference panel in that row, so damped closures appear faded.
# Output: output/vorticity_<lat>N_nx<nx>.png
ENV["GKSwstype"] = "100"            # GR without a display
using JLD2, Plots, FFTW, Statistics, Printf
gr()

const YEAR = 365 * 86400.0
const CLOSURES = ["none", "nemo", "geometric", "geometric_nosplit", "meke", "reynolds"]

length(ARGS) >= 2 || error("usage: plot_vorticity.jl <lat> <nx> [times in years...]")
lat = parse(Float64, ARGS[1])
nx  = parse(Int, ARGS[2])
times_requested = length(ARGS) > 2 ? parse.(Float64, ARGS[3:end]) : Float64[]
outdir = "output"
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

"Gaussian filter (width 2Δx_coarse) + spectral truncation, as in src/coarsegrain.jl."
function coarsegrain2d(Z, n, L; FGR = 2.0)
    nT = size(Z, 1)
    nT == n && return Float64.(Z)
    nT > n || error("cannot coarse-grain a $(nT)² field to a finer $(n)² grid")
    Zh = rfft(Float64.(Z))
    kr = 2π .* (0:nT÷2) ./ L
    l  = 2π .* fftfreq(nT, nT) ./ L
    K2 = kr .^ 2 .+ (l .^ 2)'
    Zh .*= exp.(-K2 .* (FGR * L / n)^2 ./ 24)
    lidx = vcat(1:n÷2, (nT - n÷2 + 1):nT)
    B = Zh[1:n÷2+1, lidx] .* (n^2 / nT^2)
    B[end, :] .= 0
    B[:, n÷2 + 1] .= 0
    return irfft(B, n)
end

# --- load ----------------------------------------------------------------------
columns = Tuple{String, Dict{String, Any}}[]
truthfile = find_truth(outdir, tag)
truthfile === nothing || push!(columns, ("truth (coarse-grained)", load(truthfile)))
for cl in CLOSURES
    file = joinpath(outdir, @sprintf("%s_nx%d_%s_seed1.jld2", tag, nx, cl))
    isfile(file) ? push!(columns, (cl, load(file))) : @warn "missing $file, column skipped"
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

file = joinpath(outdir, @sprintf("vorticity_%s_nx%d.png", tag, nx))
savefig(fig, file)
@info "Saved $file"
