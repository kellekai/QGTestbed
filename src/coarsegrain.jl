# ----------------------------------------------------------------------------
# Coarse-graining of truth runs and a-priori (offline) evaluation
#
# Coarse-graining operator: Gaussian filter of width Δ̄ = FGR·Δx_coarse on
# the truth grid, followed by spectral truncation to the coarse grid
# (FGR = 2 as in Perezhogin & Glazunov 2023's QG experiments).
#
# Subgrid PV forcing (what the closure should add to the coarse tendency):
#     S = J(ψ̄, q̄) − \overline{J(ψ, q)},   J(ψ, q) = ∇·(u q)
# Terms linear in the eddy fields (imposed U, β, background PV gradient,
# drag) commute with the filter and contribute nothing.
# Done on the CPU; this is an offline step.
# ----------------------------------------------------------------------------

function truncate_spectral(Ah_hi, m_hi, m_lo)
    A = Array(Ah_hi)
    nx_hi, ny_hi = m_hi.grid.nx, m_hi.grid.ny
    nx_lo, ny_lo = m_lo.grid.nx, m_lo.grid.ny
    nkr_lo, nl_hi = nx_lo ÷ 2 + 1, size(A, 2)
    lidx = vcat(1:ny_lo÷2, (nl_hi - ny_lo÷2 + 1):nl_hi)
    B = A[1:nkr_lo, lidx, :] .* (nx_lo * ny_lo / (nx_hi * ny_hi))
    B[nkr_lo, :, :] .= 0          # remove Nyquist modes of the coarse grid
    B[:, ny_lo÷2 + 1, :] .= 0
    return device_array(m_lo.grid.device)(B)
end

"Filter (Gaussian, width FGR·Δx_coarse) and truncate a truth spectral field."
function coarsegrain(Ah_hi, m_hi, m_lo; FGR = 2.0)
    G = gaussian_transfer(FGR * m_lo.Δx, m_hi)
    return truncate_spectral(G .* Ah_hi, m_hi, m_lo)
end

coarsegrain_pv(prob_hi, m_hi, m_lo; FGR = 2.0) = coarsegrain(prob_hi.sol, m_hi, m_lo; FGR)

"""
    subgrid_forcing(prob_hi, m_hi, m_lo; FGR = 2.0)

Return `(q̄h, Sh)`: coarse-grained PV and the exact subgrid PV forcing on the
coarse grid (both spectral, coarse grid, (nkr, nl, 2)).
"""
function subgrid_forcing(prob_hi, m_hi, m_lo; FGR = 2.0)
    fh = layerfields(prob_hi.sol, m_hi)
    Jh_hi = divh(fh.u .* fh.q, fh.v .* fh.q, m_hi)          # spectral ∇·(u q), truth grid
    Jbar = coarsegrain(Jh_hi, m_hi, m_lo; FGR)
    qbarh = coarsegrain(prob_hi.sol, m_hi, m_lo; FGR)
    fl = layerfields(qbarh, m_lo)
    Jres = divh(fl.u .* fl.q, fl.v .* fl.q, m_lo)
    Sh = (Jres .- Jbar) .* m_lo.mask
    return qbarh, Sh
end

"""
    apriori_skill(closure, prob_hi, m_hi, m_lo; FGR = 2.0)

Offline skill of `closure` on one truth snapshot: per-layer correlation and
relative RMS error between modelled and true subgrid forcing, plus the
energy transfer of each. For prognostic closures (GEOMETRIC, MEKE) this uses
the closure's *current* eddy-energy field, so evaluate them online, or after
spinning up the eddy energy on coarse-grained truth fields.
"""
function apriori_skill(closure, prob_hi, m_hi, m_lo; FGR = 2.0)
    qbarh, Sh = subgrid_forcing(prob_hi, m_hi, m_lo; FGR)
    f = layerfields(qbarh, m_lo)
    Fh = similar(Sh)
    pv_tendency!(Fh, closure, f, m_lo)
    S, F = tophys(Sh, m_lo), tophys(Fh, m_lo)
    corr = [cor(vec(Array(view(S, :, :, j))), vec(Array(view(F, :, :, j)))) for j in 1:2]
    rmse = [sqrt(mean((view(F, :, :, j) .- view(S, :, :, j)).^2) / mean(view(S, :, :, j).^2)) for j in 1:2]
    return (; corr, rel_rmse = rmse,
              energy_transfer_true = closure_energy_rate(Sh, f, m_lo),
              energy_transfer_model = closure_energy_rate(Fh, f, m_lo))
end
