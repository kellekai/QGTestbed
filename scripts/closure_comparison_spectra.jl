# Spectra and spectral energy transfers for the closure comparison (closure_comparison.sh), the
# diagnostics of Perezhogin & Glazunov (2023), Fig. 10, for every closure. Model environment:
#
#   julia --project=. scripts/closure_comparison_spectra.jl <outdir> [nx_coarse=64] [t_min_years=0]
#
# t_min_years: only snapshots at t ≥ t_min are used, so runs with different spin-ups (e.g. truth runs
# of the sweep, spin-up 5, linked into <outdir>) are averaged over the same period.
#
# Same definitions as perezhogin_validation_spectra.jl; every member found in <outdir>
# (truth, none, nemo, geometric, geometric_nosplit, meke, reynolds), averaged over all snapshots
# and seeds:
#   E(k)      KE spectrum, thickness-weighted over the layers, Σ_j (H_j/H) ½|k ψ̂_j|²
#   T_KE(k)   resolved KE transfer   Σ_j (H_j/H) Re[ψ̂_j* Ĵ(ψ_j, ζ_j)]
#   T_APE(k)  resolved APE transfer  (H₁F₁/H) Re[(ψ̂₁ − ψ̂₂)* Ĵ(ψ₁, ψ₂)]
#   P(k)      energy transfer of the closure, −Σ_j (H_j/H) Re[ψ̂_j* Ŝ_j], S = the closure's full PV
#             tendency on the snapshot (GM, viscosity and backscatter together; the exponential
#             filter of FilteredRK4 is not part of it). geometric/meke: recomputed with the
#             eddy-energy field stored with each snapshot.
# Truth, coarse-grained to nx_coarse: exact subgrid transfer, SGS Gauss and SGS Sharp (Ross et al. 2023
# operators 2 and 1). All densities per unit angular wavenumber, shells of width 2π/L, sum × dk =
# domain mean. Output: <outdir>/spectra.jld2 (read by plot_closure_comparison.jl).
using QGTestbed
using QGTestbed: build_problem, layerfields, divh, truncate_spectral, gaussian_transfer,
                 pv_tendency!, has_prognostic_state, set_closure_state!
using GeophysicalFlows, FourierFlows, JLD2, Printf, Statistics
using GeophysicalFlows: MultiLayerQG

outdir = length(ARGS) >= 1 ? ARGS[1] : "closure-comparison"
nx_c   = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 64
tmin   = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) * 365 * 86400.0 : 0.0
isdir(outdir) || error("folder $outdir does not exist")

const PATTERN = r"^(\d+)N_nx(\d+)_(.+)_seed(\d+)\.jld2$"
const ORDER = ["truth", "none", "nemo", "nemo_eorca025", "geometric", "geometric_nosplit", "geometric_eorca025", "meke", "reynolds"]

# --- spectral helpers (as in perezhogin_validation_spectra.jl) ---------------------------
function halfplane_weights(m)
    nkr, nl = size(m.grid.Krsq)
    w = fill(2.0, nkr, nl); w[1, :] .= 1; iseven(m.grid.nx) && (w[end, :] .= 1)
    return w
end
cross(ah, bh, m, w) = w .* real.(conj.(Array(ah)) .* Array(bh)) ./ (m.grid.nx * m.grid.ny)^2
function shells(D, m)
    K = sqrt.(Array(m.grid.Krsq)); dk = 2π / m.grid.Lx; nb = m.grid.nx ÷ 2
    s = zeros(nb + 1)
    for idx in eachindex(D)
        b = round(Int, K[idx] / dk)
        b <= nb && (s[b + 1] += D[idx])
    end
    return s ./ dk
end
wavenumbers(m) = collect((0:m.grid.nx ÷ 2) .* (2π / m.grid.Lx))

function resolved_spectra(f, m, w)
    H, Hj = m.Htot, m.H
    mask = Array(m.mask)
    E = zeros(size(m.grid.Krsq)); TKE = zero(E)
    for j in 1:2
        ψh, ζh = view(f.ψh, :, :, j), view(f.ζh, :, :, j)
        u, v, ζ = view(f.u, :, :, j), view(f.v, :, :, j), view(f.ζ, :, :, j)
        Jh = mask .* Array(divh(u .* ζ, v .* ζ, m))
        E   .+= Hj[j] / H .* (-0.5) .* cross(ψh, ζh, m, w)
        TKE .+= Hj[j] / H .* cross(ψh, Jh, m, w)
    end
    u1, v1, ψ2 = view(f.u, :, :, 1), view(f.v, :, :, 1), view(f.ψ, :, :, 2)
    J12 = mask .* Array(divh(u1 .* ψ2, v1 .* ψ2, m))
    dψh = view(f.ψh, :, :, 1) .- view(f.ψh, :, :, 2)
    TAPE = Hj[1] * m.F[1] / H .* cross(dψh, J12, m, w)
    return shells(E, m), shells(TKE, m), shells(TAPE, m)
end

function forcing_transfer(f, Sh, m, w)
    P = zeros(size(m.grid.Krsq))
    for j in 1:2
        P .-= m.H[j] / m.Htot .* cross(view(f.ψh, :, :, j), view(Sh, :, :, j), m, w)
    end
    return shells(P, m)
end

function sgs_forcing(sol_hi, m_hi, m_lo, G)
    fh = layerfields(sol_hi, m_hi)
    Jh = divh(fh.u .* fh.q, fh.v .* fh.q, m_hi)
    Jbar = truncate_spectral(G .* Jh, m_hi, m_lo)
    qbarh = truncate_spectral(G .* sol_hi, m_hi, m_lo)
    fl = layerfields(qbarh, m_lo)
    Jres = divh(fl.u .* fl.q, fl.v .* fl.q, m_lo)
    return fl, (Jres .- Jbar) .* m_lo.mask
end

function sharp_transfer(m_hi, m_lo)
    κ = sqrt.(m_hi.grid.Krsq) .* m_lo.Δx
    κc = 0.65π
    return @. ifelse(κ < κc, 1.0, exp(-23.6 * (κ - κc)^4))
end

"The closure of a run file, with the scalar settings stored in it (published values unless changed)."
function closure_from_file(member, d)
    member in ("truth", "none") && return NoClosure()
    cl = make_closure(Symbol(member))
    for (k, v) in d["closure_settings"]
        s = Symbol(k)
        (hasfield(typeof(cl), s) && !(s in (:ahm, :aei, :C_R_last))) || continue
        T = fieldtype(typeof(cl), s)
        v isa T || (v = try convert(T, v) catch; continue end)
        setfield!(cl, s, v)
    end
    return cl
end

# --- loop over runs ----------------------------------------------------------------------
files = Dict{String, Vector{String}}()
for f in sort(readdir(outdir))
    mt = match(PATTERN, f)
    mt === nothing && continue
    push!(get!(files, mt[3], String[]), joinpath(outdir, f))
end
haskey(files, "truth") || error("no truth run in $outdir")

results = Dict{String, Any}()
meta = Dict{String, Any}()
for member in [ORDER; sort([k for k in keys(files) if !(k in ORDER)])]
    haskey(files, member) || continue
    acc = Dict{String, Vector{Float64}}()
    add!(key, v) = haskey(acc, key) ? (acc[key] .+= v) : (acc[key] = copy(v))
    nsnap, nx, kgrid, kcoarse, settings, nfiles = 0, 0, Float64[], Float64[], Dict{String, Any}(), 0
    for file in files[member]
        d = load(file)
        pp, nx = d["params"], d["nx"]
        isempty(d["snapshots"]) && (@warn "$file has no PV snapshots, skipped"; continue)
        get(d, "stopped_early", false) && (@warn "$file stopped early (NaN), skipped"; continue)
        cl = closure_from_file(member, d)
        afrac = haskey(d, "aliased_fraction") ? d["aliased_fraction"] : 1/3
        prob, m, cl = build_problem(pp, nx; closure = cl, aliased_fraction = afrac)
        states = get(d, "closure_snapshots", Matrix{Float64}[])
        with_P = !(cl isa NoClosure)
        if with_P && has_prognostic_state(cl) && length(states) != length(d["snapshots"])
            @warn "$file: no eddy-energy field stored with the snapshots (run made before this script); closure transfer skipped"
            with_P = false
        end
        w = halfplane_weights(m)
        kgrid = wavenumbers(m)
        if member == "truth"
            _, m_lo, _ = build_problem(pp, nx_c; aliased_fraction = afrac)
            global w_lo = halfplane_weights(m_lo)
            kcoarse = wavenumbers(m_lo)
            global G_gauss = gaussian_transfer(2 * m_lo.Δx, m)
            global G_sharp = sharp_transfer(m, m_lo)
            global m_lo_ref = m_lo
        end
        meta["lat"], meta["Ld"], meta["L"], meta["nx_coarse"] = pp.lat, pp.Ld, pp.L, nx_c
        meta["aliased_fraction"] = afrac
        settings = Dict{String, Any}(d["closure_settings"]); settings["stepper"] = get(d, "stepper", "?")
        Fh = similar(m.cplx3)
        times = get(d, "snaptimes", fill(Inf, length(d["snapshots"])))
        used = 0
        for (i, q) in enumerate(d["snapshots"])
            times[i] >= tmin || continue
            used += 1
            MultiLayerQG.set_q!(prob, q)
            f = layerfields(prob.sol, m)
            E, TKE, TAPE = resolved_spectra(f, m, w)
            add!("E", E); add!("T_KE", TKE); add!("T_APE", TAPE)
            if with_P
                has_prognostic_state(cl) && set_closure_state!(cl, states[i])
                pv_tendency!(Fh, cl, f, m)
                add!("P", forcing_transfer(f, Fh, m, w))
            elseif member == "truth"
                fl, S = sgs_forcing(prob.sol, m, m_lo_ref, G_gauss)
                add!("P_gauss", forcing_transfer(fl, S, m_lo_ref, w_lo))
                fl, S = sgs_forcing(prob.sol, m, m_lo_ref, G_sharp)
                add!("P_sharp", forcing_transfer(fl, S, m_lo_ref, w_lo))
            end
            nsnap += 1
        end
        nfiles += 1
        @printf("%-44s %4d of %4d snapshots used (t ≥ %.1f yr), aliased_fraction %.3f, stepper %s\n", basename(file),
                used, length(d["snapshots"]), tmin / (365 * 86400.0), afrac, get(d, "stepper", "?"))
    end
    nsnap == 0 && continue
    out = Dict{String, Any}(k => v ./ nsnap for (k, v) in acc)
    out["k"], out["nx"], out["nsnap"], out["nseeds"], out["settings"] = kgrid, nx, nsnap, nfiles, settings
    member == "truth" && (out["k_coarse"] = kcoarse)
    results[member] = out
end

# --- summary: domain integrals and KE relative to the truth ------------------------------
integral(r, key) = haskey(r, key) ? sum(r[key]) * (r["k"][2] - r["k"][1]) : NaN
KEt = integral(results["truth"], "E")
@printf("\n%-20s %6s %12s %8s %12s %12s %12s\n", "member", "grid", "KE m²/s²", "KE/truth", "ΣT_KE dk", "ΣT_APE dk", "ΣP dk")
for member in [m for m in [ORDER; sort(collect(keys(results)))] if haskey(results, m)] |> unique
    r = results[member]
    @printf("%-20s %4d² %12.4e %8.3f %+12.3e %+12.3e %+12.3e\n", member, r["nx"], integral(r, "E"),
            integral(r, "E") / KEt, integral(r, "T_KE"), integral(r, "T_APE"), integral(r, "P"))
end
let t = results["truth"], dk = t["k_coarse"][2] - t["k_coarse"][1]
    @printf("%-20s %6s %12s %8s %12s %12s %+12.3e\n", "SGS Gauss (truth)", "", "", "", "", "", sum(t["P_gauss"]) * dk)
    @printf("%-20s %6s %12s %8s %12s %12s %+12.3e\n", "SGS Sharp (truth)", "", "", "", "", "", sum(t["P_sharp"]) * dk)
end

file = joinpath(outdir, "spectra.jld2")
jldsave(file; results, meta)
println("saved $file")
