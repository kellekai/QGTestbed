# Closure comparison in the layout of Perezhogin & Glazunov (2023), Fig. 10, from
# <outdir>/spectra.jld2 (written by closure_comparison_spectra.jl). Plotting environment:
#
#   julia --project=plotting scripts/plot_closure_comparison.jl closure-comparison [closures]
#
# closures: optional comma-separated subset to draw, e.g. none,nemo,reynolds (default: all found).
# (a) KE spectrum, (b) resolved KE transfer, (c) resolved APE transfer, (d) energy transfer of
# each closure (GM + viscosity + backscatter) and the exact subgrid transfer of the coarse-grained
# truth (SGS Gauss, SGS Sharp). Truth is drawn in black, no closure in gray.
# Output: <outdir>/closure_comparison.png and the curves as CSV (spectra_<member>.csv, sgs_truth.csv).
ENV["GKSwstype"] = "100"
using JLD2, Plots, Printf
gr()
include(joinpath(@__DIR__, "plotting_utils.jl"))     # closure_label

outdir = length(ARGS) >= 1 ? ARGS[1] : "closure-comparison"
d = load(joinpath(outdir, "spectra.jld2"))
R, meta = d["results"], d["meta"]
haskey(R, "truth") || error("no truth spectra in $outdir/spectra.jld2")

nxc = meta["nx_coarse"]
kLd = 1 / meta["Ld"]
kmax = π * nxc / meta["L"]
xl = (1e-5, 1.1 * kmax)
ntruth = R["truth"]["nx"]

STYLE = Dict(
    "none"              => (color = :gray50,     ls = :solid,   lw = 2.0),
    "nemo"              => (color = :darkorange, ls = :solid,   lw = 2.0),
    "nemo_eorca025"     => (color = :firebrick,  ls = :dash,    lw = 2.0),
    "geometric"         => (color = :seagreen,   ls = :dash,    lw = 2.0),
    "geometric_nosplit" => (color = :olive,      ls = :dot,     lw = 2.0),
    "geometric_eorca025" => (color = :darkgreen, ls = :dashdotdot, lw = 2.0),
    "meke"              => (color = :royalblue,  ls = :dashdot, lw = 2.0),
    "reynolds"          => (color = :purple,     ls = :dash,    lw = 2.5),
    "truth"             => (color = :black,      ls = :solid,   lw = 1.2),
)
ALL = ["none", "nemo", "nemo_eorca025", "geometric", "geometric_nosplit", "geometric_eorca025", "meke", "reynolds", "truth"]
wanted = length(ARGS) >= 2 ? [strip.(split(ARGS[2], ",")); "truth"] : ALL
ORDER = [m for m in ALL if haskey(R, m) && m in wanted]
curve_label(m) = m == "truth" ? "$(ntruth)² (truth)" : "$(nxc)², " * closure_label(m)

sel(k) = (k .> 0) .& (k .<= xl[2])
common = (xscale = :log10, xlims = xl, xlabel = "wavenumber, m⁻¹", framestyle = :box,
          guidefontsize = 10, tickfontsize = 8, titlefontsize = 11, grid = false,
          left_margin = 4Plots.mm, bottom_margin = 3Plots.mm)

function transfer_panel(key, title; legend = false)
    p = plot(; common..., title, ylabel = "×10⁻⁶ m³/s³", legend)
    hline!(p, [0.0]; color = :lightgray, ls = :dot, label = "")
    for mbr in ORDER
        r = R[mbr]; s = STYLE[mbr]; i = sel(r["k"])
        plot!(p, r["k"][i], r[key][i] .* 1e6; color = s.color, ls = s.ls, lw = s.lw, label = curve_label(mbr))
    end
    vline!(p, [kLd]; color = :steelblue, alpha = 0.6, label = "")
    return p
end

# (a) KE spectrum
pa = plot(; common..., yscale = :log10, title = "(a) KE spectrum", ylabel = "E(k), m³/s²",
          legend = :bottomleft, legendfontsize = 7)
emax = 0.0
for mbr in ORDER
    r = R[mbr]; i = sel(r["k"]) .& (r["E"] .> 0)
    global emax = max(emax, maximum(r["E"][i]))
end
for mbr in ORDER
    r = R[mbr]; s = STYLE[mbr]; i = sel(r["k"]) .& (r["E"] .> emax / 1e3)
    plot!(pa, r["k"][i], r["E"][i]; color = s.color, ls = s.ls, lw = s.lw, label = curve_label(mbr))
end
ylims!(pa, (emax / 150, emax * 2))
let k0 = 8e-5, k1 = 1.5e-4, e0 = emax * 0.8
    plot!(pa, [k0, k1], [e0, e0 * (k1 / k0)^-3]; color = :black, lw = 1.2, label = "")
    annotate!(pa, 1.25e-4, e0 * 0.9, text("k⁻³", 9, :left))
end
vline!(pa, [kLd]; color = :steelblue, alpha = 0.6, label = "")
annotate!(pa, kLd * 1.05, emax / 60, text(@sprintf("deformation\nradius %.0f km", meta["Ld"] / 1e3), 8, :steelblue, :left))

pb = transfer_panel("T_KE", "(b) Resolved KE transfer")
pc = transfer_panel("T_APE", "(c) Resolved APE transfer")

# (d) closure transfer vs exact subgrid transfer
pd = plot(; common..., title = "(d) Subgrid KE+APE transfer", ylabel = "×10⁻⁶ m³/s³",
          legend = :bottomleft, legendfontsize = 7)
hline!(pd, [0.0]; color = :lightgray, ls = :dot, label = "")
for mbr in ORDER
    r = R[mbr]; haskey(r, "P") || continue
    s = STYLE[mbr]; i = sel(r["k"])
    plot!(pd, r["k"][i], r["P"][i] .* 1e6; color = s.color, ls = s.ls, lw = s.lw, label = curve_label(mbr))
end
let t = R["truth"], kc = t["k_coarse"], i = (kc .> 0)
    plot!(pd, kc[i], t["P_gauss"][i] .* 1e6; color = :black, ls = :dashdot, lw = 1.3, label = "SGS Gauss ($(ntruth)² → $(nxc)²)")
    plot!(pd, kc[i], t["P_sharp"][i] .* 1e6; color = :black, ls = :dash,    lw = 1.3, label = "SGS Sharp ($(ntruth)² → $(nxc)²)")
end
vline!(pd, [kLd]; color = :steelblue, alpha = 0.6, label = "")

# Regime as in the proposal's Ld/Δx map: < 1 eddies unresolved, 1–2 gray zone, > 2 resolved
Ld_dx = meta["Ld"] / (meta["L"] / nxc)
regime = Ld_dx < 1 ? "unresolved" : Ld_dx <= 2 ? "gray zone" : "resolved"
fig = plot(pa, pb, pc, pd; layout = (2, 2), size = (1150, 880),
           plot_title = @sprintf("Closures at %.0f°N, %d² (Ld/Δx = %.2f, %s), truth %d², Ld = %.1f km",
                                 meta["lat"], nxc, Ld_dx, regime, ntruth, meta["Ld"] / 1e3),
           plot_titlefontsize = 12)

# CSV copies of the curves
function write_csv(path, cols)
    open(path, "w") do io
        println(io, join(first.(cols), ","))
        for i in eachindex(last(cols[1]))
            println(io, join([@sprintf("%.8e", c[i]) for c in last.(cols)], ","))
        end
    end
end
for mbr in ORDER
    r = R[mbr]
    cols = ["k" => r["k"], "E" => r["E"], "T_KE" => r["T_KE"], "T_APE" => r["T_APE"]]
    haskey(r, "P") && push!(cols, "P" => r["P"])
    write_csv(joinpath(outdir, "spectra_$(mbr).csv"), cols)
end
write_csv(joinpath(outdir, "sgs_truth.csv"),
          ["k_coarse" => R["truth"]["k_coarse"], "P_gauss" => R["truth"]["P_gauss"], "P_sharp" => R["truth"]["P_sharp"]])

name = length(ARGS) >= 2 ? "closure_comparison_" * replace(ARGS[2], "," => "-") * ".png" : "closure_comparison.png"
file = joinpath(outdir, name)
savefig(fig, file)
println("saved $file")
for mbr in ORDER
    @printf("  %-28s %2d seed(s), %4d snapshots\n", curve_label(mbr), R[mbr]["nseeds"], R[mbr]["nsnap"])
end
