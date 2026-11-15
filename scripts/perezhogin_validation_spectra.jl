# Spectra and spectral energy transfers for the Perezhogin & Glazunov (2023) validation
# (their Fig. 8), from the PV snapshots of the runs in <outdir>. Model environment:
#
#   julia --project=. scripts/perezhogin_validation_spectra.jl perezhogin-validation [nx_coarse=64]
#
# For every run (truth, none, reynolds; averaged over all snapshots and seeds):
#   E(k)      KE spectrum, thickness-weighted over the layers, Σ_j (H_j/H) ½|k ψ̂_j|²
#   T_KE(k)   resolved KE transfer   Σ_j (H_j/H) Re[ψ̂_j* Ĵ(ψ_j, ζ_j)]
#   T_APE(k)  resolved APE transfer  (H₁F₁/H) Re[(ψ̂₁ − ψ̂₂)* Ĵ(ψ₁, ψ₂)]
#   P(k)      subgrid transfer of total energy by the closure, −Σ_j (H_j/H) Re[ψ̂_j* Ŝ_j]
#             (reynolds only; S = the closure's PV tendency on the snapshot)
# and for the truth, coarse-grained to nx_coarse, the transfer of the exact subgrid forcing
#   S = J(ψ̄, q̄) − \overline{J(ψ, q)} with two operators of Ross et al. (2023):
#   P_gauss   operator 2: Gaussian exp(−k²(2Δx)²/24) + spectral truncation ("SGS Gauss")
#   P_sharp   operator 1: pyqg's filter exp(−23.6(κ* − 0.65π)⁴) + truncation ("SGS Sharp")
# All are densities per unit wavenumber in shells of width 2π/L, normalized so that the sum
# times dk is the domain mean (m³/s² for E, m³/s³ for the transfers); k is the angular
# wavenumber (rad/m), as in pyqg. Ĵ includes the model's 2/3 dealiasing mask.
# Output: <outdir>/spectra.jld2 (read by plot_perezhogin_validation.jl).
using QGTestbed
using QGTestbed: build_problem, layerfields, divh, truncate_spectral, gaussian_transfer,
                 pv_tendency!, ReynoldsBackscatter
using GeophysicalFlows, FourierFlows, JLD2, Printf, Statistics
using GeophysicalFlows: MultiLayerQG

outdir = length(ARGS) >= 1 ? ARGS[1] : "perezhogin-validation"
nx_c   = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 64
isdir(outdir) || error("folder $outdir does not exist")

const PATTERN = r"^(\d+)N_nx(\d+)_(truth|none|reynolds)_seed(\d+)\.jld2$"

# --- spectral helpers ------------------------------------------------------------------
function halfplane_weights(m)
    nkr, nl = size(m.grid.Krsq)
    w = fill(2.0, nkr, nl); w[1, :] .= 1; iseven(m.grid.nx) && (w[end, :] .= 1)
    return w
end
"Domain-mean density of Re(â* b̂) per mode (2D spectral fields)."
cross(ah, bh, m, w) = w .* real.(conj.(Array(ah)) .* Array(bh)) ./ (m.grid.nx * m.grid.ny)^2
"Sum a per-mode density into isotropic shells of width dk = 2π/L; returns density per dk."
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

"E, T_KE, T_APE of one state on its own grid."
function resolved_spectra(f, m, w)
    H, Hj = m.Htot, m.H
    mask = Array(m.mask)
    E = zeros(size(m.grid.Krsq)); TKE = zero(E)
    for j in 1:2
        ψh, ζh = view(f.ψh, :, :, j), view(f.ζh, :, :, j)
        u, v, ζ = view(f.u, :, :, j), view(f.v, :, :, j), view(f.ζ, :, :, j)
        Jh = mask .* Array(divh(u .* ζ, v .* ζ, m))                 # Ĵ(ψ_j, ζ_j) = ∇·(u_j ζ_j)
        E   .+= Hj[j] / H .* (-0.5) .* cross(ψh, ζh, m, w)          # ½|kψ̂|² = −½Re(ψ̂* ζ̂)
        TKE .+= Hj[j] / H .* cross(ψh, Jh, m, w)
    end
    u1, v1, ψ2 = view(f.u, :, :, 1), view(f.v, :, :, 1), view(f.ψ, :, :, 2)
    J12 = mask .* Array(divh(u1 .* ψ2, v1 .* ψ2, m))                # Ĵ(ψ₁, ψ₂)
    dψh = view(f.ψh, :, :, 1) .- view(f.ψh, :, :, 2)
    TAPE = Hj[1] * m.F[1] / H .* cross(dψh, J12, m, w)
    return shells(E, m), shells(TKE, m), shells(TAPE, m)
end

"Energy transfer −Σ (H_j/H) Re[ψ̂_j* Ŝ_j] of a PV tendency S (spectral, (nkr, nl, 2))."
function forcing_transfer(f, Sh, m, w)
    P = zeros(size(m.grid.Krsq))
    for j in 1:2
        P .-= m.H[j] / m.Htot .* cross(view(f.ψh, :, :, j), view(Sh, :, :, j), m, w)
    end
    return shells(P, m)
end

"Exact subgrid PV forcing on the coarse grid for a filter transfer function G on the truth grid."
function sgs_forcing(sol_hi, m_hi, m_lo, G)
    fh = layerfields(sol_hi, m_hi)
    Jh = divh(fh.u .* fh.q, fh.v .* fh.q, m_hi)
    Jbar = truncate_spectral(G .* Jh, m_hi, m_lo)
    qbarh = truncate_spectral(G .* sol_hi, m_hi, m_lo)
    fl = layerfields(qbarh, m_lo)
    Jres = divh(fl.u .* fl.q, fl.v .* fl.q, m_lo)
    return fl, (Jres .- Jbar) .* m_lo.mask
end

"pyqg's scale-selective filter of the coarse grid, evaluated on the truth grid (Ross et al. operator 1)."
function sharp_transfer(m_hi, m_lo)
    κ = sqrt.(m_hi.grid.Krsq) .* m_lo.Δx
    κc = 0.65π
    return @. ifelse(κ < κc, 1.0, exp(-23.6 * (κ - κc)^4))
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
for member in ("truth", "none", "reynolds")
    haskey(files, member) || (@warn "no '$member' runs in $outdir, skipped"; continue)
    acc = Dict{String, Vector{Float64}}()
    add!(key, v) = haskey(acc, key) ? (acc[key] .+= v) : (acc[key] = copy(v))
    nsnap, nx, kgrid, kcoarse, CR = 0, 0, Float64[], Float64[], NaN
    for file in files[member]
        d = load(file)
        pp, nx = d["params"], d["nx"]
        isempty(d["snapshots"]) && (@warn "$file has no PV snapshots, skipped"; continue)
        cl = member == "reynolds" ? ReynoldsBackscatter(; C_R = Float64(d["closure_settings"]["C_R"])) : NoClosure()
        member == "reynolds" && (CR = cl.C_R)
        afrac = haskey(d, "aliased_fraction") ? d["aliased_fraction"] : 1/3   # as in the run
        prob, m, cl = build_problem(pp, nx; closure = cl, aliased_fraction = afrac)
        w = halfplane_weights(m)
        kgrid = wavenumbers(m)
        if member == "truth"
            _, m_lo, _ = build_problem(pp, nx_c; aliased_fraction = afrac)
            w_lo = halfplane_weights(m_lo)
            kcoarse = wavenumbers(m_lo)
            G_gauss = gaussian_transfer(2 * m_lo.Δx, m)
            G_sharp = sharp_transfer(m, m_lo)
        end
        meta["lat"], meta["Ld"], meta["L"], meta["nx_coarse"] = pp.lat, pp.Ld, pp.L, nx_c
        meta["aliased_fraction"] = afrac
        for q in d["snapshots"]
            MultiLayerQG.set_q!(prob, q)
            f = layerfields(prob.sol, m)
            E, TKE, TAPE = resolved_spectra(f, m, w)
            add!("E", E); add!("T_KE", TKE); add!("T_APE", TAPE)
            if member == "reynolds"
                Fh = similar(m.cplx3)
                pv_tendency!(Fh, cl, f, m)
                add!("P_sgs", forcing_transfer(f, Fh, m, w))
            elseif member == "truth"
                fl, S = sgs_forcing(prob.sol, m, m_lo, G_gauss)
                add!("P_gauss", forcing_transfer(fl, S, m_lo, w_lo))
                fl, S = sgs_forcing(prob.sol, m, m_lo, G_sharp)
                add!("P_sharp", forcing_transfer(fl, S, m_lo, w_lo))
            end
            nsnap += 1
        end
        @printf("%-40s %4d snapshots\n", basename(file), length(d["snapshots"]))
    end
    nsnap == 0 && continue
    out = Dict{String, Any}(k => v ./ nsnap for (k, v) in acc)
    out["k"], out["nx"], out["nsnap"], out["nseeds"] = kgrid, nx, nsnap, length(files[member])
    member == "truth" && (out["k_coarse"] = kcoarse)
    member == "reynolds" && (out["C_R"] = CR)
    results[member] = out
    dk = kgrid[2] - kgrid[1]
    @printf("%-8s %dx%d  KE = %.4e m²/s²  ΣT_KE·dk = %+.2e  ΣT_APE·dk = %+.2e m²/s³\n", member, nx, nx,
            sum(out["E"]) * dk, sum(out["T_KE"]) * dk, sum(out["T_APE"]) * dk)
end

file = joinpath(outdir, "spectra.jld2")
jldsave(file; results, meta)
println("saved $file")
