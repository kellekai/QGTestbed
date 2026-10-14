# ----------------------------------------------------------------------------
# Reynolds-stress backscatter (Perezhogin & Glazunov 2023, JAMES,
# arXiv:2304.06789), as used in their 2-layer QG experiments.
#
# Germano decomposition of the subfilter PV flux: Leonard + Cross + Reynolds.
# The Reynolds part (eddy–eddy) carries the KE backscatter and is modelled by
# scale similarity on the resolved small scales (their Eq. 19):
#
#   q′ = q̄ − G(q̄),  u′ = ū − G(ū)
#   σ_j^KEB = G(u′_j q′) − G(u′_j) G(q′)
#   ∂q/∂t  += −C_R ∂_j σ_j^KEB                     (their Eq. 25, layerwise)
#
# G is a Gaussian filter, transfer exp(−k²Δ̄²/24), with Δ̄ = 2Δx in the QG
# experiments (filter-to-grid ratio 2). The same C_R in both layers.
#
# In their QG runs only this term is added on top of pyqg's scale-selective
# spectral filter (here: FilteredRK4); the dissipative DMM part was used in
# the barotropic and NEMO experiments, not in QG. C_R = 7 was tuned to the
# large-scale KE spectrum at 64² (15.6 km) in the pyqg "eddy" configuration,
# which is the 49°N member of this sweep. Keep C_R fixed across latitudes to
# test the non-scale-aware behaviour, or set `dynamic_CR = true` to use
# their energy-balance estimate (Eqs. 21–22) with σ^DMM replaced by the
# resolved Leonard-type flux, see the docstring of `dynamic_CR`.
# ----------------------------------------------------------------------------

Base.@kwdef mutable struct ReynoldsBackscatter <: AbstractClosure
    C_R::Float64 = 7.0
    FGR::Float64 = 2.0               # Δ̄ / Δx
    dynamic_CR::Bool = false
    C_R_last::Float64 = NaN
end

default_stepper(::ReynoldsBackscatter) = "FilteredRK4"

function keb_flux(f, Δf, m)
    Gq, Gu, Gv = gfilter(f.q, Δf, m), gfilter(f.u, Δf, m), gfilter(f.v, Δf, m)
    qp, up, vp = f.q .- Gq, f.u .- Gu, f.v .- Gv
    Gqp, Gup, Gvp = gfilter(qp, Δf, m), gfilter(up, Δf, m), gfilter(vp, Δf, m)
    σx = gfilter(up .* qp, Δf, m) .- Gup .* Gqp
    σy = gfilter(vp .* qp, Δf, m) .- Gvp .* Gqp
    return σx, σy
end

"""
    dynamic_CR(f, σx, σy, Δf, m)

Experimental. Energy-balance estimate in the spirit of PG23 Eqs. 21–22:
choose C_R such that the total modelled flux σ = σ^L + C_R σ^KEB satisfies
⟨σ·β⟩ = 0 with β = ∇ψ − (Δ̄²/12)∇q, where σ^L = G(ū q̄) − G(ū)G(q̄) is the
computable Leonard-type flux standing in for the paper's σ^DMM. The paper
does not use the dynamic C_R in QG; compare against fixed C_R = 7 at 49°N
before relying on it.
"""
function dynamic_CR(f, σx, σy, Δf, m)
    Gq, Gu, Gv = gfilter(f.q, Δf, m), gfilter(f.u, Δf, m), gfilter(f.v, Δf, m)
    σLx = gfilter(f.u .* f.q, Δf, m) .- Gu .* Gq
    σLy = gfilter(f.v .* f.q, Δf, m) .- Gv .* Gq
    βx = ∂x(f.ψ, m) .- (Δf^2 / 12) .* ∂x(f.q, m)
    βy = ∂y(f.ψ, m) .- (Δf^2 / 12) .* ∂y(f.q, m)
    num = sum(σLx .* βx .+ σLy .* βy)
    den = sum(σx .* βx .+ σy .* βy)
    return den == 0 ? 0.0 : max(-num / den, 0.0)
end

function pv_tendency!(Fqh, c::ReynoldsBackscatter, f, m)
    Δf = c.FGR * m.Δx
    σx, σy = keb_flux(f, Δf, m)
    CR = c.dynamic_CR ? dynamic_CR(f, σx, σy, Δf, m) : c.C_R
    c.C_R_last = CR
    Fqh .= -CR .* divh(σx, σy, m)
    Fqh .*= m.mask
    return Fqh
end

"Domain-mean meridional PV flux of the closure: tendency −C_R∇·σ ⇒ flux C_R σ."
function param_pv_flux(c::ReynoldsBackscatter, f, m)
    _, σy = keb_flux(f, c.FGR * m.Δx, m)
    CR = c.dynamic_CR ? c.C_R_last : c.C_R
    return (CR * mean(view(σy, :, :, 1)), CR * mean(view(σy, :, :, 2)))
end

closure_diagnostics(c::ReynoldsBackscatter, m) = (C_R = c.C_R_last,)
