# Scale-awareness summary: every coarse run reduced to three numbers relative to
# the truth at the same latitude, plotted against Δx/Ld.
#
#   julia --project=plotting scripts/scale_awareness.jl            # all runs in output/
#   julia --project=plotting scripts/scale_awareness.jl -s 5 -d output
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
# row per grid (32², 64², 128²) and the same Δx/Ld axis in every row:
#   <dir>/scale_awareness_spectral_error.png
#   <dir>/scale_awareness_resolved_KE.png
#   <dir>/scale_awareness_PV_flux.png
ENV["GKSwstype"] = "100"
using JLD2, Plots, Printf
gr()
include(joinpath(@__DIR__, "plotting_utils.jl"))

# --- arguments -------------------------------------------------------------------
outdir, spinup = "output", 5.0
let i = 1
    while i <= length(ARGS)
        if ARGS[i] == "-d"
            global outdir = ARGS[i+1]; i += 2
        elseif ARGS[i] == "-s"
            global spinup = parse(Float64, ARGS[i+1]); i += 2
        else
            error("unknown argument $(ARGS[i])")
        end
    end
end

const CLOSURE_ORDER = ["none", "nemo", "geometric", "geometric_nosplit", "meke", "reynolds"]
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

runs = RunSummary[]
for file in sort(readdir(outdir))
    m = match(PATTERN, file)
    m === nothing && continue
    lat, nx, closure, seed = parse(Int, m[1]), parse(Int, m[2]), m[3], parse(Int, m[4])
    r = summarize(joinpath(outdir, file), lat, nx, closure, seed)
    r === nothing || push!(runs, r)
end
isempty(runs) && error("no usable run files in $outdir")

# --- truth reference per latitude (averaged over truth seeds) -----------------------
truth = Dict{Int, NamedTuple}()
for lat in unique(r.lat for r in runs)
    trs = [r for r in runs if r.lat == lat && r.closure == "truth"]
    isempty(trs) && (@warn "no truth run for $(lat)°N: latitude skipped"; continue)
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
    println(io, "lat,nx,closure,dx_over_Ld,n_seeds,spectral_error,spectral_error_min,spectral_error_max,resolved_KE_rel_error,resolved_KE_min,resolved_KE_max,PV_flux_rel_error,PV_flux_min,PV_flux_max")
    for key in sort(collect(keys(groups)))
        g = groups[key]
        r = g[1].r
        col(f) = [getfield(x.m, f) for x in g]
        s, e, q = col(:spec_err), col(:ke_ratio), col(:flux_ratio)
        @printf(io, "%d,%d,%s,%.4f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
                key[1], key[2], key[3], (r.L / r.nx) / r.Ld, length(g),
                mean(s), minimum(s), maximum(s), mean(e), minimum(e), maximum(e), mean(q), minimum(q), maximum(q))
    end
end
@info "Saved $csvfile"

# --- plots: one figure per metric, one row per grid ----------------------------------
closures = [c for c in CLOSURE_ORDER if any(k[3] == c for k in keys(groups))]
append!(closures, sort(unique(k[3] for k in keys(groups) if !(k[3] in CLOSURE_ORDER))))
grids = sort(unique(k[2] for k in keys(groups)))
dxld(k) = (groups[k][1].r.L / k[2]) / groups[k][1].r.Ld
allx = [dxld(k) for k in keys(groups)]
xl = (minimum(allx) / 1.15, maximum(allx) * 1.15)       # same x-range in every row

const METRICS = [
    (:spec_err,   "spectral error (RMS log₁₀ E/E_truth, k ≤ k_c/2)", "spectral_error"),
    (:ke_ratio,   "resolved KE / truth − 1",                          "resolved_KE"),
    (:flux_ratio, "upper-layer PV flux / truth − 1",                  "PV_flux"),
]

for (f, title, fname) in METRICS
    rows = Plots.Plot[]
    for (ri, nx) in enumerate(grids)
        p = plot(xscale = :log10, xlims = xl, legend = false, titlefontsize = 10, guidefontsize = 9,
                 tickfontsize = 8, title = "$(nx)² grid", ylabel = ri == cld(length(grids), 2) ? title : "",
                 xlabel = ri == length(grids) ? "Δx / Ld" : "")
        hline!(p, [0.0]; color = :black, lw = 0.5)
        for (ci, cl) in enumerate(closures)
            ks = sort([k for k in keys(groups) if k[3] == cl && k[2] == nx]; by = dxld)
            isempty(ks) && continue
            x = dxld.(ks)
            y = [mean(getfield(g.m, f) for g in groups[k]) for k in ks]
            lo = [y[i] - minimum(getfield(g.m, f) for g in groups[ks[i]]) for i in eachindex(ks)]
            hi = [maximum(getfield(g.m, f) for g in groups[ks[i]]) - y[i] for i in eachindex(ks)]
            plot!(p, x, y; yerror = (lo, hi), color = ci, marker = :circle, markersize = 4, lw = 1.5, label = "")
        end
        push!(rows, p)
    end
    # legend column: closures in the first row, empty space below
    leg = plot(framestyle = :none, legend = :left, legendfontsize = 9, grid = false, showaxis = false)
    for (ci, cl) in enumerate(closures)
        plot!(leg, [NaN], [NaN]; color = ci, lw = 2, marker = :circle, label = cl)
    end
    blank = plot(framestyle = :none, grid = false, showaxis = false, legend = false)
    panels = Plots.Plot[]
    for (ri, p) in enumerate(rows)
        push!(panels, p, ri == 1 ? leg : blank)
    end
    fig = plot(panels...; layout = grid(length(grids), 2, widths = [0.84, 0.16]),
               size = (1100, 330 * length(grids)), left_margin = 8Plots.mm, bottom_margin = 4Plots.mm)
    png = joinpath(outdir, "scale_awareness_$(fname).png")
    savefig(fig, png)
    @info "Saved $png"
end
