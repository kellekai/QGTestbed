# Score a tuning experiment against the truth and find the best parameter set.
#
#   julia --project=plotting scripts/tune_evaluate.jl --closure reynolds
#   julia --project=plotting scripts/tune_evaluate.jl --closure geometric --weights 1,1,0 --spinup 3
#
# Options
#   --closure NAME     closure that was tuned (reads output/tune/<closure>/ unless --dir is given)
#   --dir DIR          tuning output directory
#   --truthdir DIR     where the truth runs are (default output)
#   --spinup S         spin-up of the tuning runs excluded from the scores (default 3)
#   --truth-spinup S   spin-up excluded from the truth (default 5)
#   --weights a,b,c    cost = a·spectral error + b·|resolved KE error| + c·|PV-flux error|
#                      (default 1,1,1; same metrics as scale_awareness.jl)
#
# For every parameter set the cost is averaged over all cases (latitude:grid) and
# seeds. If the as-published run of the same closure exists in --truthdir for a
# case, its cost is shown as the reference ("published").
#
# Output in the tuning directory:
#   tune_results.csv   one row per parameter set and case, plus the mean over cases, ranked
#   tune_cost.png      cost against the parameter(s): a line for one parameter, a map for two
ENV["GKSwstype"] = "100"
using JLD2, Plots, Printf
gr()
include(joinpath(@__DIR__, "plotting_utils.jl"))

closure, dir, truthdir, spinup, truth_spinup, weights = "", "", "output", 3.0, 5.0, [1.0, 1.0, 1.0]
let i = 1
    while i <= length(ARGS)
        a, v = ARGS[i], (i < length(ARGS) ? ARGS[i+1] : error("option $(ARGS[i]) needs a value"))
        if a == "--closure";          global closure = v
        elseif a == "--dir";          global dir = v
        elseif a == "--truthdir";     global truthdir = v
        elseif a == "--spinup";       global spinup = parse(Float64, v)
        elseif a == "--truth-spinup"; global truth_spinup = parse(Float64, v)
        elseif a == "--weights";      global weights = parse.(Float64, split(v, ","))
        else error("unknown option $a")
        end
        i += 2
    end
end
isempty(closure) && isempty(dir) && error("give --closure (or --dir)")
dir = isempty(dir) ? joinpath("output", "tune", closure) : dir
length(weights) == 3 || error("--weights needs three numbers")

# --- score every tuning run ---------------------------------------------------------
truths = Dict{Int, Any}()
rows = NamedTuple[]
for file in sort(readdir(dir))
    endswith(file, ".jld2") || continue
    r = load_summary(joinpath(dir, file); spinup)
    r === nothing && (@warn "$file: incomplete, skipped"; continue)
    isempty(r.metadata) && (@warn "$file: no tuning parameters stored, skipped"; continue)
    lat = round(Int, r.lat)
    t = get!(truths, lat) do
        truth_reference(truthdir, lat; spinup = truth_spinup)
    end
    t === nothing && (@warn "no truth run for $(lat)°N in $truthdir; $file skipped"; continue)
    m = run_metrics(r, t)
    push!(rows, (; lat, nx = r.nx, params = Dict(String(k) => Float64(v) for (k, v) in r.metadata),
                   dx_Ld = (r.L / r.nx) / r.Ld, m..., cost = cost(m, weights)))
end
isempty(rows) && error("no scored tuning runs in $dir")

pnames = sort(unique(vcat([collect(keys(r.params)) for r in rows]...)))
pkey(r) = Tuple(get(r.params, n, NaN) for n in pnames)
cases = sort(unique((r.lat, r.nx) for r in rows))

# mean over seeds per (parameter set, case), then over cases per parameter set
percase = Dict{Tuple, Vector{Float64}}()       # (pkey, case) => [cost, spec, ke, flux]
for g in unique((pkey(r), (r.lat, r.nx)) for r in rows)
    sel = [r for r in rows if pkey(r) == g[1] && (r.lat, r.nx) == g[2]]
    percase[g] = [mean(getfield(r, f) for r in sel) for f in (:cost, :spec_err, :ke_ratio, :flux_ratio)]
end
psets = sort(unique(pkey(r) for r in rows))
overall = Dict(p => mean(v[1] for (g, v) in percase if g[1] == p) for p in psets)
ncases = Dict(p => count(g -> g[1] == p, keys(percase)) for p in psets)
ranked = sort(psets; by = p -> overall[p])

# reference: the as-published run of the same closure, if present for the cases
published = Dict{Tuple{Int, Int}, Float64}()
for (lat, nx) in cases
    f = joinpath(truthdir, @sprintf("%02dN_nx%d_%s_seed1.jld2", lat, nx, closure))
    isfile(f) || continue
    r = load_summary(f; spinup = truth_spinup)
    r === nothing && continue
    published[(lat, nx)] = cost(run_metrics(r, truths[lat]), weights)
end

# --- table ---------------------------------------------------------------------------
csv = joinpath(dir, "tune_results.csv")
open(csv, "w") do io
    println(io, join(vcat("rank", pnames, ["case", "dx_over_Ld", "cost", "spectral_error", "resolved_KE_rel_error", "PV_flux_rel_error"]), ","))
    for (rank, p) in enumerate(ranked)
        vals = [@sprintf("%g", x) for x in p]
        println(io, join(vcat(string(rank), vals, ["all cases", "", @sprintf("%.4f", overall[p]), "", "", ""]), ","))
        for c in cases
            haskey(percase, (p, c)) || continue
            v = percase[(p, c)]
            dxl = first(r.dx_Ld for r in rows if (r.lat, r.nx) == c)
            println(io, join(vcat(string(rank), vals, [@sprintf("%dN %d²", c...), @sprintf("%.3f", dxl),
                                  (@sprintf("%.4f", x) for x in v)...]), ","))
        end
    end
end

println("\nTuning of $closure: $(length(psets)) parameter sets, cases ", join(["$(c[1])°N/$(c[2])²" for c in cases], ", "),
        ", weights (spectrum, KE, flux) = ", weights)
!isempty(published) && println("published settings, cost per case: ",
                               join([@sprintf("%d°N/%d²: %.3f", c..., v) for (c, v) in sort(collect(published))], ", "),
                               @sprintf("  (mean %.3f)", mean(values(published))))
println("best parameter sets (mean cost over cases):")
for p in ranked[1:min(5, end)]
    @printf("  %s   cost = %.3f%s\n", join(["$n = $(@sprintf("%g", x))" for (n, x) in zip(pnames, p)], ", "),
            overall[p], ncases[p] < length(cases) ? "  (only $(ncases[p]) of $(length(cases)) cases finished)" : "")
end
println("table: $csv")

# --- plot ----------------------------------------------------------------------------
png = joinpath(dir, "tune_cost.png")
if length(pnames) == 1
    x = [p[1] for p in psets]
    order = sortperm(x)
    logx = all(>(0), x) && maximum(x) / minimum(x) > 20
    fig = plot(xlabel = pnames[1], ylabel = "cost", xscale = logx ? :log10 : :identity, legend = :best,
               title = "$closure: cost against $(pnames[1])", titlefontsize = 10, size = (800, 450))
    for (ci, c) in enumerate(cases)
        y = [get(percase, (p, c), [NaN])[1] for p in psets]
        plot!(fig, x[order], y[order]; marker = :circle, color = ci, label = "$(c[1])°N, $(c[2])²")
        haskey(published, c) && hline!(fig, [published[c]]; color = ci, linestyle = :dash, label = "published, $(c[1])°N")
    end
    plot!(fig, x[order], [overall[p] for p in psets][order]; color = :black, lw = 2, marker = :diamond, label = "mean over cases")
    savefig(fig, png); @info "Saved $png"
elseif length(pnames) == 2
    xs, ys = sort(unique(p[1] for p in psets)), sort(unique(p[2] for p in psets))
    Z = [get(overall, (x, y), NaN) for y in ys, x in xs]
    best = ranked[1]
    fig = heatmap(string.(xs), string.(ys), Z; xlabel = pnames[1], ylabel = pnames[2], c = :viridis,
                  colorbar_title = "mean cost", title = "$closure: mean cost over cases (★ = best)",
                  titlefontsize = 10, size = (800, 550))
    scatter!(fig, [findfirst(==(best[1]), xs) - 0.5], [findfirst(==(best[2]), ys) - 0.5];
             marker = :star5, markersize = 12, color = :white, label = "")
    savefig(fig, png); @info "Saved $png"
else
    @info "More than two parameters: no plot, see $csv"
end
