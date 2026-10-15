# Energy time series and time-mean upper-layer KE spectra for one or more runs, overlaid.
#
#   julia --project=plotting scripts/plot_summary.jl output/50N_nx64_none_seed1.jld2
#   julia --project=plotting scripts/plot_summary.jl output/40N_nx64_*_seed1.jld2 output/40N_nx256_truth_seed1.jld2 -o output/summary_40N_nx64.png
#
# Left: upper-layer KE (solid) and eddy PE (dashed) against time, with the
# default 5-year spin-up marked; use it to check that each run has equilibrated.
# Right: time-mean upper-layer KE spectrum E(k), averaged over the vorticity
# snapshots after the spin-up and normalized so that ∫E(k)dk = domain-mean KE.
# The normalization is independent of the grid, so truth and coarse runs are
# directly comparable. Dotted vertical lines (labelled 'grid cutoff') mark the
# largest wavenumber kept by the dealiasing (2/3 rule) of each coarse grid.
#
# Options (after the files): `-o name.png` output file, `-s 5` spin-up in years.
ENV["GKSwstype"] = "100"            # GR without a display
using JLD2, Plots, Printf
gr()
include(joinpath(@__DIR__, "plotting_utils.jl"))     # YEAR, ke_spectrum_from_vorticity, mean_spectrum, dealias_cutoff

# --- arguments -------------------------------------------------------------------
files, outfile, spinup = String[], nothing, 5.0
let i = 1
    while i <= length(ARGS)
        if ARGS[i] == "-o"
            global outfile = ARGS[i+1]; i += 2
        elseif ARGS[i] == "-s"
            global spinup = parse(Float64, ARGS[i+1]); i += 2
        else
            push!(files, ARGS[i]); i += 1
        end
    end
end
isempty(files) && error("give one or more .jld2 files")
outfile = something(outfile, joinpath("output", "summary_" * splitext(basename(files[1]))[1] * ".png"))

# --- plot ------------------------------------------------------------------------
pE = plot(xlabel = "years", ylabel = "energy [m²/s²]", yscale = :log10, legend = :bottomright,
          legendfontsize = 8, title = "KE upper (solid), PE (dashed)", titlefontsize = 9)
vline!(pE, [spinup]; color = :gray, linestyle = :dot, label = @sprintf("spin-up (%.0f yr)", spinup))
pS = plot(xlabel = "k [rad/m]", ylabel = "upper-layer KE spectrum [m³/s²]", xscale = :log10, yscale = :log10,
          legend = :bottomleft, legendfontsize = 8, titlefontsize = 9,
          title = @sprintf("time mean after %.0f yr", spinup))

allspec = Float64[]
lats = Dict{Float64, Float64}()           # latitude => Ld (m)
cutoffs = Dict{Int, Float64}()          # nx => dealiasing cutoff
for (i, file) in enumerate(files)
    d = load(file)
    lats[d["params"].lat] = d["params"].Ld
    label = @sprintf("%s %d²", closure_label(d["closure"]), d["nx"])
    ts = d["timeseries"]
    t = ts["t"] ./ YEAR
    for (key, ls, lab) in (("KE1", :solid, label), ("PE", :dash, ""))
        y = ts[key]; ok = y .> 0
        plot!(pE, t[ok], y[ok]; color = i, linestyle = ls, label = lab)
    end

    if !haskey(d, "vorticity_upper")
        @warn "$file has no vorticity snapshots; spectrum skipped"
        continue
    end
    L = d["params"].L
    res = mean_spectrum(d["vorticity_upper"], d["vorticity_times"], L, spinup)
    if res === nothing
        @warn "$file has no snapshots after $(spinup) yr; spectrum skipped"
        continue
    end
    k, S = res
    ok = (k .> 0) .& (S .> 0)
    plot!(pS, k[ok], S[ok]; color = i, label = label)
    append!(allspec, S[ok])
    d["closure"] == "truth" || (cutoffs[d["nx"]] = dealias_cutoff(d["nx"], L))
end

for (nx, kc) in sort(collect(cutoffs))
    vline!(pS, [kc]; color = :gray, linestyle = :dot, label = "grid cutoff $(nx)²")
end
if !isempty(allspec)                    # keep decayed runs from stretching the axis
    top = maximum(allspec)
    ylims!(pS, (top * 1e-8, top * 5))
end

title = join([@sprintf("%.0f°N (Ld = %.1f km)", φ, Ld / 1e3) for (φ, Ld) in sort(collect(lats))], ", ")
fig = plot(pE, pS; layout = (1, 2), size = (1100, 450), left_margin = 5Plots.mm, bottom_margin = 5Plots.mm,
           plot_title = title, plot_titlefontsize = 12)
savefig(fig, outfile)
@info "Saved $outfile"
