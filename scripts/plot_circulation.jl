# Mean flow and eddies, in the style of Jansen et al. (2019), Figure 3, for the truth
# and every closure at one latitude and grid.
#
#   julia --project=plotting scripts/plot_circulation.jl 40 64          # snapshot = last saved time
#   julia --project=plotting scripts/plot_circulation.jl 40 64 12.5     # snapshot at 12.5 years
#   julia --project=plotting scripts/plot_circulation.jl 40 64 12.5 -s 5
#
# Contours: sea surface height. In QG the surface pressure is ρ₀ f₀ ψ₁, so
#   η = f₀ ψ₁ / g,
# with ψ₁ the upper-layer streamfunction including the imposed flow (ψ₁ = ψ₁' − U₁ y).
# Black contours positive, white negative.
# Shading: magnitude of the vertically integrated transport |H₁u₁ + H₂u₂| (m²/s),
# imposed flow included, on a logarithmic colour scale. Runs saved before the
# lower-layer streamfunction was stored (no `psi_lower`) fall back to the
# upper-layer speed |u₁| (m/s); the figure title says which one is shown.
#
# Row (a): time mean after the spin-up (default 5 yr).
# Row (b): snapshot of the deviations from that time mean (SSH and transport).
# Columns: truth (coarse-grained to the same grid), then each closure.
# Contour intervals and colour ranges are shared within a row, set from the truth.
#
# In this doubly periodic, homogeneous testbed the time mean is mostly the imposed
# flow (straight zonal SSH contours, uniform transport); departures from it are
# eddy-driven zonal jets and sampling noise. Row (b) is the informative one.
#
# Output: output/circulation_<lat>N_nx<nx>.png
ENV["GKSwstype"] = "100"
using JLD2, Plots, Printf
gr()
include(joinpath(@__DIR__, "plotting_utils.jl"))

const G = 9.81
const USE_LOWER = Ref(true)          # set to false when not every run has psi_lower
const CLOSURES = ["none", "nemo", "geometric", "geometric_nosplit", "meke", "reynolds"]

# --- arguments ---------------------------------------------------------------------
pos = String[]; spinup = 5.0
let i = 1
    while i <= length(ARGS)
        if ARGS[i] == "-s"
            global spinup = parse(Float64, ARGS[i+1]); i += 2
        else
            push!(pos, ARGS[i]); i += 1
        end
    end
end
length(pos) >= 2 || error("usage: plot_circulation.jl <lat> <nx> [snapshot time in years] [-s spinup]")
lat, nx = parse(Float64, pos[1]), parse(Int, pos[2])
tsnap = length(pos) >= 3 ? parse(Float64, pos[3]) : nothing
outdir = "output"
tag = @sprintf("%02.0fN", lat)

# --- load runs ------------------------------------------------------------------------
columns = Tuple{String, String}[]
truthfiles = filter(f -> startswith(f, tag * "_nx") && endswith(f, "_truth_seed1.jld2"), readdir(outdir))
isempty(truthfiles) || push!(columns, ("truth (coarse-grained)", joinpath(outdir, first(truthfiles))))
for cl in CLOSURES
    f = joinpath(outdir, @sprintf("%s_nx%d_%s_seed1.jld2", tag, nx, cl))
    isfile(f) ? push!(columns, (cl, f)) : @warn "missing $f, column skipped"
end
isempty(columns) && error("no runs for $tag at $(nx)² in $outdir")

"Time-mean and snapshot-deviation fields (SSH, transport components) on the coarse grid."
function circulation(file)
    jldopen(file, "r") do f
        pp = f["params"]
        L, f0 = pp.L, pp.f₀
        Z, T = f["vorticity_upper"], f["vorticity_times"]
        have2 = USE_LOWER[] && haskey(f, "psi_lower") && !isempty(f["psi_lower"])
        P2 = have2 ? f["psi_lower"] : nothing
        sel = findall(>(spinup * YEAR), T)
        isempty(sel) && error("$file: no snapshots after the spin-up ($spinup yr)")
        js = tsnap === nothing ? lastindex(T) : argmin(abs.(T ./ YEAR .- tsnap))
        n = size(Z[1], 1)
        y = collect(range(-L / 2, step = L / n, length = n))

        # eddy fields of one snapshot → (ψ1', and either transport or upper-layer velocity)
        function fields(j)
            ψ1 = psi_from_vorticity(Z[j], L)
            u1, v1 = velocities(ψ1, L)
            if have2
                u2, v2 = velocities(P2[j], L)
                return ψ1, pp.H1 .* u1 .+ pp.H2 .* u2, pp.H1 .* v1 .+ pp.H2 .* v2
            else
                return ψ1, u1, v1
            end
        end
        # time mean of the eddy parts
        ψm = zeros(n, n); am = zeros(n, n); bm = zeros(n, n)
        for j in sel
            ψ, a, b = fields(j)
            ψm .+= ψ; am .+= a; bm .+= b
        end
        ψm ./= length(sel); am ./= length(sel); bm ./= length(sel)
        ψs, as, bs = fields(js)

        # coarse-grain (truth) to the plotted grid
        cg(A) = coarsegrain2d(A, nx, L)
        ψm, am, bm, ψs, as, bs = cg.((ψm, am, bm, ψs, as, bs))
        yc = collect(range(-L / 2, step = L / nx, length = nx))

        # add the imposed flow to the mean: ψ1 = ψ1' − U1 y ; transport x-component + H1 U1 (or U1)
        Ubg = have2 ? pp.H1 * pp.U[1] + pp.H2 * pp.U[2] : pp.U[1]
        ψmean = ψm .- pp.U[1] .* yc'
        mean_ssh = f0 .* ψmean ./ G
        mean_mag = sqrt.((am .+ Ubg) .^ 2 .+ bm .^ 2)
        dev_ssh = f0 .* (ψs .- ψm) ./ G
        dev_mag = sqrt.((as .- am) .^ 2 .+ (bs .- bm) .^ 2)
        return (; mean_ssh, mean_mag, dev_ssh, dev_mag, have2, tshow = T[js] / YEAR, L, Ld = pp.Ld)
    end
end

data = [(name, circulation(f)) for (name, f) in columns]
have2 = all(d.have2 for (_, d) in data)
any(d.have2 for (_, d) in data) && !have2 &&
    @warn "some runs lack psi_lower: all columns show upper-layer speed for consistency"
if !have2 && any(d.have2 for (_, d) in data)
    # recompute with the upper-layer fallback for all runs, so the columns are comparable
    USE_LOWER[] = false
    data = [(name, circulation(f)) for (name, f) in columns]
end
qty = have2 ? "|vertically integrated transport| [m²/s]" : "|upper-layer velocity| [m/s]"

L = data[1][2].L
xkm = collect(range(-L / 2, step = L / nx, length = nx)) ./ 1e3

"Round a positive number to 1, 2 or 5 × 10ⁿ."
function nice(x)
    e = floor(log10(x)); m = x / 10^e
    return (m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10) * 10^e
end

function panel(mag, ssh, crange, cint, title, ylabel, xlabel, showcb)
    p = heatmap(xkm, xkm, permutedims(log10.(max.(mag, 10.0^crange[1])));
                c = :viridis, clims = crange, aspect_ratio = :equal, colorbar = false,
                xlims = extrema(xkm), ylims = extrema(xkm), framestyle = :box,
                title = title, titlefontsize = 9, tickfontsize = 6, guidefontsize = 7,
                xlabel = xlabel, ylabel = ylabel)
    lo, hi = extrema(ssh)
    pos_levels = collect(cint:cint:max(hi, cint))
    neg_levels = collect(-cint:-cint:min(lo, -cint))
    isempty(pos_levels) || contour!(p, xkm, xkm, permutedims(ssh); levels = pos_levels, color = :black, lw = 0.8, colorbar = false)
    isempty(neg_levels) || contour!(p, xkm, xkm, permutedims(ssh); levels = neg_levels, color = :white, lw = 0.8, colorbar = false)
    return p
end

panels = Plots.Plot[]
ncol = length(data)
for (row, (magkey, sshkey, label)) in enumerate(((:mean_mag, :mean_ssh, "(a) time mean"),
                                                   (:dev_mag, :dev_ssh, "(b) snapshot − time mean")))
    ref = data[1][2]
    refmag = getfield(ref, magkey)
    crange = (log10(max(quantile(vec(refmag), 0.02), 1e-12)), log10(maximum(refmag)))
    crange[2] > crange[1] || (crange = (crange[1], crange[1] + 1))
    sshref = getfield(ref, sshkey)
    cint = nice((maximum(sshref) - minimum(sshref)) / 10 + eps())
    for (c, (name, d)) in enumerate(data)
        t = row == 1 ? (c == 1 ? "$name\n$label" : name) :
                       (c == 1 ? @sprintf("%s, %.1f yr", label, d.tshow) : @sprintf("%.1f yr", d.tshow))
        push!(panels, panel(getfield(d, magkey), getfield(d, sshkey), crange, cint, t,
                            c == 1 ? "y [km]" : "", row == 2 ? "x [km]" : "", false))
    end
    # colour bar of this row as a narrow strip
    lv = collect(range(crange[1], crange[2], length = 201))
    push!(panels, heatmap([0.0, 1.0], lv, repeat(lv, 1, 2); c = :viridis, clims = crange, colorbar = false,
                          xticks = false, ymirror = true, tickfontsize = 6, framestyle = :box,
                          title = @sprintf("log₁₀\nSSH ci %.2g m", cint), titlefontsize = 7,
                          right_margin = 8Plots.mm))
end

wpanel = 0.95 / ncol
fig = plot(panels...; layout = grid(2, ncol + 1, widths = [fill(wpanel, ncol); 1 - ncol * wpanel]),
           size = (270 * ncol + 100, 620),
           plot_title = @sprintf("%.0f°N, %d² (Δx/Ld = %.2f): SSH contours (black > 0, white < 0), shading %s",
                                 lat, nx, (L / nx) / data[1][2].Ld, qty),
           plot_titlefontsize = 10)
file = joinpath(outdir, @sprintf("circulation_%s_nx%d.png", tag, nx))
savefig(fig, file)
@info "Saved $file"
