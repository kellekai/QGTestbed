# Scale-awareness summary: every coarse run reduced to three numbers relative to
# the truth at the same latitude, plotted against the resolution parameter Ld/Δx
# (grid points per deformation radius; Hallberg 2013). Shading: Ld/Δx < 1 unresolved,
# 1–2 gray zone, > 2 resolved.
#
#   julia --project=plotting scripts/scale_awareness.jl            # all runs in output/
#   julia --project=plotting scripts/scale_awareness.jl -s 5 -d output
#   julia --project=plotting scripts/scale_awareness.jl -d output/tuned -t output   # tuned sweep
#
# -d DIR   folder with the coarse runs (and output folder of the table and figures)
# -t DIR   folder with the truth runs (default: same as -d); if the -d folder has no
#          no-closure runs ('none'), they are also taken from here as the reference
# -s S     spin-up in years excluded from the time means (default 5)
# -g N,…   grids to plot, e.g. -g 64 or -g 32,64 (default: all); the figure names then
#          end in _nx64 (or _nx32-64), the CSV always contains all grids
#
#   julia --project=plotting scripts/scale_awareness.jl -g 64     # only the 64² grid
#
# Metrics (all time means after the spin-up, default 5 yr):
#
#   spectral error   RMS of log10(E_model / E_truth) over the well-resolved scales,
#                    0 < k ≤ k_c/2, with k_c the dealiasing cutoff of the coarse grid.
#                    0 = perfect; 1 = off by a factor 10 on average. Spectra are the
#                    upper-layer KE spectra from the vorticity snapshots (grid independent).
#                    A run whose spectrum lies more than 6 decades below the truth is
#                    counted as 6 (eddies gone).
#   resolved KE      ∫_0^{k_c} E_model dk / ∫_0^{k_c} E_truth dk − 1
#                    (the truth's energy at scales the coarse grid can represent).
#   PV flux          ⟨vq⟩_resolved + ⟨vq⟩_param in the upper layer, divided by the
#                    truth's ⟨vq⟩, minus 1. For GM-type closures the parametrized
#                    part is included (it is diagnosed only, see README).
#
# With several seeds the mean is plotted, with the range as an error bar. The
# truth spectrum and flux are averaged over all available truth seeds.
#
# Output: <dir>/scale_awareness.csv and three figures, one per metric, with one
# row per grid (32², 64², 128²), a linear Ld/Δx axis fitted to each grid's range (bottom)
# and the latitude of each point on a second x-axis on top of each row:
#   <dir>/scale_awareness_spectral_error.png
#   <dir>/scale_awareness_resolved_KE.png
#   <dir>/scale_awareness_PV_flux.png
ENV["GKSwstype"] = "100"
using JLD2, Plots, Printf
gr()
include(joinpath(@__DIR__, "plotting_utils.jl"))

# --- arguments -------------------------------------------------------------------
outdir, spinup, truthdir, gridsel = "output", 5.0, "", Int[]
let i = 1
    while i <= length(ARGS)
        if ARGS[i] == "-d"
            global outdir = ARGS[i+1]; i += 2
        elseif ARGS[i] == "-s"
            global spinup = parse(Float64, ARGS[i+1]); i += 2
        elseif ARGS[i] == "-t"
            global truthdir = ARGS[i+1]; i += 2
        elseif ARGS[i] == "-g"
            global gridsel = parse.(Int, split(ARGS[i+1], ",")); i += 2
        else
            error("unknown argument $(ARGS[i])")
        end
    end
end

const CLOSURE_ORDER = ["none", "nemo", "geometric", "geometric_nosplit", "meke", "reynolds"]
# Legend names come from LABELS in plotting_utils.jl.
const MAX_DECADES = 6.0
const PATTERN = r"^(\d+)N_nx(\d+)_(.+)_seed(\d+)\.jld2$"

# --- read every run once (only what is needed) -------------------------------------
struct RunSummary
    lat::Int
    nx::Int
    closure::String
    seed::Int
    L::Float64
    Ld::Float64
    k::Vector{Float64}
    S::Vector{Float64}
    vq1::Float64          # resolved + parametrized, upper layer, time mean
end

function summarize(path, lat, nx, closure, seed)
    jldopen(path, "r") do f
        pp = f["params"]
        haskey(f, "vorticity_upper") || (@warn "$path: no vorticity snapshots, skipped"; return nothing)
        res = mean_spectrum(f["vorticity_upper"], f["vorticity_times"], pp.L, spinup)
        res === nothing && (@warn "$path: no snapshots after spin-up, skipped"; return nothing)
        ts = f["timeseries"]
        late = ts["t"] .> spinup * YEAR
        haskey(ts, "vq_resolved_1") || (@warn "$path: no PV-flux time series, skipped"; return nothing)
        vq = mean(ts["vq_resolved_1"][late]) + mean(ts["vq_param_1"][late])
        return RunSummary(lat, nx, closure, seed, pp.L, pp.Ld, res[1], res[2], vq)
    end
end

truthdir = isempty(truthdir) ? outdir : truthdir
runs = RunSummary[]
function read_runs!(runs, dir, keep)
    for file in sort(readdir(dir))
        m = match(PATTERN, file)
        m === nothing && continue
        lat, nx, closure, seed = parse(Int, m[1]), parse(Int, m[2]), m[3], parse(Int, m[4])
        keep(closure) || continue
        r = summarize(joinpath(dir, file), lat, nx, closure, seed)
        r === nothing || push!(runs, r)
    end
end
if truthdir == outdir
    read_runs!(runs, outdir, c -> true)
else
    read_runs!(runs, outdir, c -> c != "truth")
    # truth from the truth folder; the no-closure reference too if the run folder has none
    none_here = any(r.closure == "none" for r in runs)
    read_runs!(runs, truthdir, c -> c == "truth" || (c == "none" && !none_here))
    none_here || any(r.closure == "none" for r in runs) &&
        @info "no-closure runs ('none') taken from $truthdir as reference"
end
isempty(runs) && error("no usable run files in $outdir")

# --- truth reference per latitude (averaged over truth seeds) -----------------------
truth = Dict{Int, NamedTuple}()
for lat in unique(r.lat for r in runs)
    trs = [r for r in runs if r.lat == lat && r.closure == "truth"]
    isempty(trs) && (@warn "no truth run for $(lat)°N in $truthdir: latitude skipped"; continue)
    nmin = minimum(length(r.S) for r in trs)
    truth[lat] = (k = trs[1].k[1:nmin], S = mean(r.S[1:nmin] for r in trs), vq1 = mean(r.vq1 for r in trs))
end

# --- metrics per run ------------------------------------------------------------------
function metrics(r::RunSummary, t)
    kc = dealias_cutoff(r.nx, r.L)
    n = min(length(r.S), length(t.S))
    half = [i for i in 2:n if r.k[i] <= kc / 2 && t.S[i] > 0]
    full = [i for i in 2:n if r.k[i] <= kc]
    floor_ = t.S[half] .* 10.0^(-MAX_DECADES)
    spec_err = sqrt(mean(log10.(max.(r.S[half], floor_) ./ t.S[half]) .^ 2))
    dk = r.k[2] - r.k[1]
    ke_ratio = sum(r.S[full]) * dk / (sum(t.S[full]) * dk) - 1
    flux_ratio = r.vq1 / t.vq1 - 1
    return (; spec_err, ke_ratio, flux_ratio)
end

rows = []
for r in runs
    r.closure == "truth" && continue
    haskey(truth, r.lat) || continue
    push!(rows, (r = r, m = metrics(r, truth[r.lat])))
end
isempty(rows) && error("no coarse runs with a matching truth run")

# group over seeds
groups = Dict{Tuple{Int, Int, String}, Vector{Any}}()
for row in rows
    push!(get!(groups, (row.r.lat, row.r.nx, row.r.closure), Any[]), row)
end

csvfile = joinpath(outdir, "scale_awareness.csv")
open(csvfile, "w") do io
    println(io, "lat,nx,closure,Ld_over_dx,dx_over_Ld,n_seeds,spectral_error,spectral_error_min,spectral_error_max,resolved_KE_rel_error,resolved_KE_min,resolved_KE_max,PV_flux_rel_error,PV_flux_min,PV_flux_max")
    for key in sort(collect(keys(groups)))
        g = groups[key]
        r = g[1].r
        col(f) = [getfield(x.m, f) for x in g]
        s, e, q = col(:spec_err), col(:ke_ratio), col(:flux_ratio)
        @printf(io, "%d,%d,%s,%.4f,%.4f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
                key[1], key[2], key[3], r.Ld / (r.L / r.nx), (r.L / r.nx) / r.Ld, length(g),
                mean(s), minimum(s), maximum(s), mean(e), minimum(e), maximum(e), mean(q), minimum(q), maximum(q))
    end
end
@info "Saved $csvfile"

# --- plots: one figure per metric, one row per grid ----------------------------------
closures = [c for c in CLOSURE_ORDER if any(k[3] == c for k in keys(groups))]
append!(closures, sort(unique(k[3] for k in keys(groups) if !(k[3] in CLOSURE_ORDER))))
grids = sort(unique(k[2] for k in keys(groups)))
if !isempty(gridsel)
    missing_grids = setdiff(gridsel, grids)
    isempty(missing_grids) || error("no runs on grid(s) $(missing_grids); available: $(grids)")
    grids = sort(gridsel)
end
suffix = isempty(gridsel) ? "" : "_nx" * join(grids, "-")
"Resolution parameter Ld/Δx (grid points per deformation radius) of a group."
ldx(k) = groups[k][1].r.Ld / (groups[k][1].r.L / k[2])
"x-range of one row: the Ld/Δx values of that grid plus a 5 % margin on each side."
function row_xlims(nx)
    xs = [ldx(k) for k in keys(groups) if k[2] == nx]
    pad = 0.05 * max(maximum(xs) - minimum(xs), 0.1 * maximum(xs))
    return (max(0.0, minimum(xs) - pad), maximum(xs) + pad)
end

const REGIMES = [(0.0, 1.0, :red), (1.0, 2.0, :orange), (2.0, Inf, :green)]
const REGIME_LABELS = ["Ld/Δx < 1: unresolved", "Ld/Δx 1–2: gray zone", "Ld/Δx > 2: resolved"]

const METRICS = [
    (:spec_err,   "spectral error (RMS log₁₀ E/E_truth, k ≤ k_c/2)", "spectral_error"),
    (:ke_ratio,   "resolved KE / truth − 1",                          "resolved_KE"),
    (:flux_ratio, "upper-layer PV flux / truth − 1",                  "PV_flux"),
]

for (f, title, fname) in METRICS
    rows = Plots.Plot[]
    for (ri, nx) in enumerate(grids)
        xl = row_xlims(nx)
        first_row = ri == 1
        p = plot(xlims = xl, guidefontsize = 9, tickfontsize = 8,
                 legend = first_row ? :topright : false, legendfontsize = 10,
                 background_color_legend = Plots.RGBA(1, 1, 1, 0.85), foreground_color_legend = :gray,
                 ylabel = ri == cld(length(grids), 2) ? title : "",
                 xlabel = ri == length(grids) ? "Ld / Δx  (grid points per deformation radius)" : "",
                 annotations = (:topleft, text("  $(nx)² grid", 10, :left)))
        # resolution regimes (Hallberg 2013): Ld/Δx < 1 unresolved, 1–2 gray zone, > 2 resolved
        for ((a, b, c), lab) in zip(REGIMES, REGIME_LABELS)
            lo_, hi_ = max(a, xl[1]), min(b, xl[2])
            lo_ < hi_ && vspan!(p, [lo_, hi_]; color = c, alpha = 0.12, lw = 0, label = first_row ? lab : "")
        end
        hline!(p, [0.0]; color = :black, lw = 0.5, label = "")
        # latitudes of the points in this row: thin guide lines and a second x-axis on top
        lk = sort(unique([(round(ldx(k); digits = 6), k[1]) for k in keys(groups) if k[2] == nx]))
        vline!(p, first.(lk); color = :gray, lw = 0.5, alpha = 0.4, linestyle = :dot, label = "")
        top = twiny(p)
        plot!(top; xlims = xl, xticks = (first.(lk), ["$(l)°N" for l in last.(lk)]),
              xlabel = ri == 1 ? "latitude" : "", yticks = false, grid = false, legend = false,
              tickfontsize = 7, guidefontsize = 9)
        for (ci, cl) in enumerate(closures)
            ks = sort([k for k in keys(groups) if k[3] == cl && k[2] == nx]; by = ldx)
            isempty(ks) && continue
            x = ldx.(ks)
            y = [mean(getfield(g.m, f) for g in groups[k]) for k in ks]
            lo = [y[i] - minimum(getfield(g.m, f) for g in groups[ks[i]]) for i in eachindex(ks)]
            hi = [maximum(getfield(g.m, f) for g in groups[ks[i]]) - y[i] for i in eachindex(ks)]
            plot!(p, x, y; yerror = (lo, hi), color = ci, marker = :circle, markersize = 4, lw = 1.5,
                  label = first_row ? closure_label(cl) : "")
        end
        push!(rows, p)
    end
    fig = plot(rows...; layout = (length(grids), 1),
               size = (1000, length(grids) == 1 ? 520 : 330 * length(grids)),
               left_margin = 8Plots.mm, bottom_margin = 4Plots.mm, right_margin = 4Plots.mm)
    png = joinpath(outdir, "scale_awareness_$(fname)$(suffix).png")
    savefig(fig, png)
    @info "Saved $png"
end
