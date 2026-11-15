# Energy budget of the MEKE closure (Jansen et al. 2019) from one or more `meke` runs.
#
#   julia --project=plotting scripts/meke_budget.jl output/40N_nx64_meke_seed1.jld2 output/60N_nx64_meke_seed1.jld2
#   julia --project=plotting scripts/meke_budget.jl output/*N_nx64_meke_seed1.jld2 -t output -o output/meke_budget_nx64.png
#
# Needs runs made with the current code (time series meke_src_GM, …; older runs: rerun them).
# Prints, per run, the time means after the spin-up of
#   e (MEKE), √(2e), K_GM, L_mix, R and the truth's κ_eff = −⟨v₁q₁⟩/(F₁ΔU) for comparison,
#   the four budget terms of ∂e/∂t (m² s⁻³) and their share of the total source,
#   the drag rate c_d/H √(e + |u_bot|² + U_bg²) as a damping time, and the part due to U_bg.
# The figure shows the terms against time (sources > 0, sinks < 0), one panel per run.
#
# Options (after the files): -t DIR truth folder for κ_eff (default: folder of the first file),
#                            -s 5 spin-up in years, -o name.png figure file.
ENV["GKSwstype"] = "100"
using JLD2, Plots, Printf, Statistics
gr()
include(joinpath(@__DIR__, "plotting_utils.jl"))     # YEAR, closure_label, load_summary

files, truthdir, spinup, outfile = String[], "", 5.0, ""
let i = 1
    while i <= length(ARGS)
        a = ARGS[i]
        if a in ("-t", "-s", "-o")
            i < length(ARGS) || error("option $a needs a value")
            a == "-t" && (global truthdir = ARGS[i+1])
            a == "-s" && (global spinup = parse(Float64, ARGS[i+1]))
            a == "-o" && (global outfile = ARGS[i+1])
            i += 2
        else
            push!(files, a); i += 1
        end
    end
end
isempty(files) && error("usage: meke_budget.jl <meke run files...> [-t DIR] [-s years] [-o figure.png]")
truthdir = isempty(truthdir) ? dirname(files[1]) : truthdir
outfile = isempty(outfile) ? joinpath(dirname(files[1]), "meke_budget.png") : outfile

const TERMS = [("meke_src_GM", "GM PE release", 1), ("meke_src_biharm", "biharmonic Smagorinsky", 1),
               ("meke_sink_bs", "backscatter", -1), ("meke_sink_drag", "drag", -1)]

function truth_kappa(dir, lat)
    isdir(dir) || return NaN
    pat = Regex("^$(round(Int, lat))N_nx\\d+_truth_seed\\d+\\.jld2\$")
    fs = filter(f -> occursin(pat, f), readdir(dir))
    isempty(fs) && return NaN
    vals = Float64[]
    for f in fs
        s = load_summary(joinpath(dir, f); spinup)
        s === nothing && continue
        p = load(joinpath(dir, f), "params")
        push!(vals, -s.vq1 / (p.F1 * p.ΔU))
    end
    return isempty(vals) ? NaN : mean(vals)
end

panels = Plots.Plot[]
for file in files
    d = load(file)
    d["closure"] == "meke" || (@warn "$file is not a meke run, skipped"; continue)
    ts = d["timeseries"]
    haskey(ts, "meke_src_GM") ||
        (@warn "$file has no MEKE budget time series (made with older code): rerun it"; continue)
    p = d["params"]; nx = d["nx"]
    t = ts["t"] ./ YEAR
    late = t .> spinup
    any(late) || (@warn "$file has no output after $(spinup) yr, skipped"; continue)
    tm(k) = mean(filter(isfinite, ts[k][late]))

    e, K, R, Lmix = tm("MEKE"), tm("kappa_gm"), tm("R"), tm("meke_Lmix")
    src = tm("meke_src_GM") + tm("meke_src_biharm")
    drag_rate, drag_Ubg = tm("meke_drag_rate"), tm("meke_drag_rate_Ubg")
    κ_eff = truth_kappa(truthdir, p.lat)

    println("="^78)
    @printf("%s   (%.0f°N, %d², Ld/Δx = %.2f)\n", basename(file), p.lat, nx, p.Ld / (p.L / nx))
    @printf("  e = %.3e m²/s²   √(2e) = %.4f m/s   L_mix = %.1f km (Δx = %.1f km)   R = %.3f\n",
            e, sqrt(2e), Lmix / 1e3, p.L / nx / 1e3, R)
    @printf("  K_GM = %.1f m²/s   truth κ_eff = %.1f m²/s   ratio = %.3f\n", K, κ_eff, K / κ_eff)
    println("  budget of ∂e/∂t, time mean after $(spinup) yr (m² s⁻³; share of total source):")
    for (k, name, sgn) in TERMS
        v = tm(k)
        @printf("    %-6s %-24s %+11.3e   %6.1f %%\n", sgn > 0 ? "source" : "sink", name, sgn * v, 100 * v / src)
    end
    @printf("    %-31s %+11.3e   (≈ d⟨e⟩/dt, should be ≈ 0)\n", "net", tm("meke_net"))
    @printf("  drag damping time 1/rate = %.1f d, of which U_bg alone gives %.1f d (U_bg share of rate %.0f %%)\n",
            1 / drag_rate / 86400, 1 / drag_Ubg / 86400, 100 * drag_Ubg / drag_rate)

    pl = plot(title = @sprintf("%.0f°N, %d² (Ld/Δx = %.2f)", p.lat, nx, p.Ld / (p.L / nx)), titlefontsize = 9,
              xlabel = "years", ylabel = "m² s⁻³", legend = :outerright, legendfontsize = 7, tickfontsize = 7)
    for (i, (k, name, sgn)) in enumerate(TERMS)
        plot!(pl, t, sgn .* ts[k]; color = i, label = (sgn > 0 ? "+ " : "− ") * name)
    end
    hline!(pl, [0.0]; color = :black, lw = 0.5, label = "")
    vline!(pl, [spinup]; color = :gray, linestyle = :dot, label = "spin-up")
    push!(panels, pl)
end
println("="^78)

if !isempty(panels)
    fig = plot(panels...; layout = (length(panels), 1), size = (900, 260 * length(panels)),
               left_margin = 6Plots.mm, bottom_margin = 3Plots.mm)
    mkpath(dirname(abspath(outfile)))
    savefig(fig, outfile)
    @info "Saved $outfile"
end
