# Figure in the layout of Perezhogin & Glazunov (2023), Fig. 8, from
# <outdir>/spectra.jld2 (written by perezhogin_validation_spectra.jl). Plotting environment:
#
#   julia --project=plotting scripts/plot_perezhogin_validation.jl perezhogin-validation
#
# (a) KE spectrum, (b) resolved KE transfer, (c) resolved APE transfer, (d) subgrid transfer of
# total energy: the Reynolds closure and the exact subgrid forcing of the coarse-grained truth
# (SGS Gauss, SGS Sharp). Runs: 64² without closure, 64² Reynolds, 256² truth (the paper's
# JansenHeld run is not part of this validation).
# Output: <outdir>/perezhogin_validation.png and the curves as CSV (spectra_<member>.csv, sgs_truth.csv)
ENV["GKSwstype"] = "100"
using JLD2, Plots, Printf
gr()

outdir = length(ARGS) >= 1 ? ARGS[1] : "perezhogin-validation"
d = load(joinpath(outdir, "spectra.jld2"))
R, meta = d["results"], d["meta"]
haskey(R, "truth") || error("no truth spectra in $outdir/spectra.jld2")

nxc = meta["nx_coarse"]
kLd = 1 / meta["Ld"]                         # the paper marks k = 1/r_d
kmax = π * nxc / meta["L"]                   # Nyquist wavenumber of the coarse grid
xl = (1e-5, 1.1 * kmax)
ntruth = R["truth"]["nx"]

STYLE = Dict(
    "none"     => (label = "$(nxc)²",            color = :gray50,  ls = :solid, lw = 2.0),
    "reynolds" => (label = "$(nxc)², Reynolds",  color = :purple,  ls = :dash,  lw = 2.5),
    "truth"    => (label = "$(ntruth)²",         color = :black,   ls = :solid, lw = 1.0),
)
ORDER = [m for m in ("none", "reynolds", "truth") if haskey(R, m)]

sel(k) = (k .> 0) .& (k .<= xl[2])
common = (xscale = :log10, xlims = xl, xlabel = "wavenumber, m⁻¹", framestyle = :box,
          guidefontsize = 10, tickfontsize = 8, titlefontsize = 11, grid = false,
          left_margin = 4Plots.mm, bottom_margin = 3Plots.mm)

function transfer_panel(key, title; ylab = "×10⁻⁶ m³/s³", legend = false)
    p = plot(; common..., title, ylabel = ylab, legend)
    hline!(p, [0.0]; color = :lightgray, ls = :dot, label = "")
    for mbr in ORDER
        r = R[mbr]; s = STYLE[mbr]; i = sel(r["k"])
        plot!(p, r["k"][i], r[key][i] .* 1e6; color = s.color, ls = s.ls, lw = s.lw, label = s.label)
    end
    vline!(p, [kLd]; color = :steelblue, alpha = 0.6, label = "")
    return p
end

# (a) KE spectrum
pa = plot(; common..., yscale = :log10, title = "(a) KE spectrum", ylabel = "E(k), m³/s²",
          legend = :bottomleft, legendfontsize = 8)
emax = 0.0
for mbr in ORDER
    r = R[mbr]; s = STYLE[mbr]; i = sel(r["k"]) .& (r["E"] .> 0)
    plot!(pa, r["k"][i], r["E"][i]; color = s.color, ls = s.ls, lw = s.lw, label = s.label)
    global emax = max(emax, maximum(r["E"][i]))
end
ylims!(pa, (emax / 150, emax * 2))
let k0 = 8e-5, k1 = 1.5e-4, e0 = emax * 0.8     # k⁻³ reference slope
    plot!(pa, [k0, k1], [e0, e0 * (k1 / k0)^-3]; color = :black, lw = 1.2, label = "")
    annotate!(pa, 1.25e-4, e0 * 0.9, text("k⁻³", 9, :left))
end
vline!(pa, [kLd]; color = :steelblue, alpha = 0.6, label = "")
annotate!(pa, kLd * 1.05, emax / 60, text(@sprintf("deformation\nradius %.0f km", meta["Ld"] / 1e3), 8, :steelblue, :left))

pb = transfer_panel("T_KE", "(b) Resolved KE transfer")
pc = transfer_panel("T_APE", "(c) Resolved APE transfer")

# (d) subgrid transfer of total energy
pd = plot(; common..., title = "(d) Subgrid KE+APE transfer", ylabel = "×10⁻⁶ m³/s³",
          legend = :bottomleft, legendfontsize = 8)
hline!(pd, [0.0]; color = :lightgray, ls = :dot, label = "")
if haskey(R, "reynolds") && haskey(R["reynolds"], "P_sgs")
    r = R["reynolds"]; s = STYLE["reynolds"]; i = sel(r["k"])
    plot!(pd, r["k"][i], r["P_sgs"][i] .* 1e6; color = s.color, ls = s.ls, lw = s.lw,
          label = @sprintf("%s (C_R = %g)", s.label, r["C_R"]))
end
let t = R["truth"], kc = t["k_coarse"], i = (kc .> 0)
    plot!(pd, kc[i], t["P_gauss"][i] .* 1e6; color = :black, ls = :dashdot, lw = 1.3, label = "SGS Gauss ($(ntruth)² → $(nxc)²)")
    plot!(pd, kc[i], t["P_sharp"][i] .* 1e6; color = :black, ls = :dash,    lw = 1.3, label = "SGS Sharp ($(ntruth)² → $(nxc)²)")
end
vline!(pd, [kLd]; color = :steelblue, alpha = 0.6, label = "")

nseeds = join(["$(STYLE[m].label): $(R[m]["nseeds"]) seed(s), $(R[m]["nsnap"]) snapshots" for m in ORDER], "; ")
fig = plot(pa, pb, pc, pd; layout = (2, 2), size = (1050, 820),
           plot_title = @sprintf("Perezhogin & Glazunov (2023) Fig. 10 setup: %.0f°N, Ld = %.1f km, %s",
                                 meta["lat"], meta["Ld"] / 1e3,
                                 get(meta, "aliased_fraction", 1/3) == 0 ? "no 2/3 dealiasing (pyqg-like)" : "2/3 dealiasing"),
           plot_titlefontsize = 12)
# CSV copies of the curves (for overlays, e.g. pyqg-validation/analyze_pyqg.py --testbed <outdir>)
function write_csv(path, cols)          # cols: vector of "name" => values
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
    haskey(r, "P_sgs") && push!(cols, "P_sgs" => r["P_sgs"])
    write_csv(joinpath(outdir, "spectra_$(mbr).csv"), cols)
end
write_csv(joinpath(outdir, "sgs_truth.csv"),
          ["k_coarse" => R["truth"]["k_coarse"], "P_gauss" => R["truth"]["P_gauss"], "P_sharp" => R["truth"]["P_sharp"]])

file = joinpath(outdir, "perezhogin_validation.png")
savefig(fig, file)
println("saved $file  ($nseeds)")
