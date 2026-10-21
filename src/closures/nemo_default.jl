# ----------------------------------------------------------------------------
# NEMO default (not scale aware)
#
# Momentum: bilaplacian viscosity with NEMO 4 grid-scaled coefficient
#     ahm = (1/12) U_v Δx³          (ldfdyn, nn_ahm_ijk_t = 20, rn_Uv)
# Eddy-induced advection: GM with a fixed, grid-scaled coefficient
#     aei = (1/2) U_e L_e            (ldftra, nn_aei_ijk_t = 20, rn_Ue, rn_Le)
#
# Defaults: U_v = 0.1 m/s (NEMO 4 namelist_ref rn_Uv; Mak et al. 2023 used 0.5 in their
# eddy-permitting runs, see geometric.jl), U_e = 0.02 m/s,
# L_e = Δx. Set them to the values in your NEMO 4.0.7 namelist
# (namdyn_ldf / namtra_eiv). Set `gm = false` if your 25 km
# configuration runs without GM.
# ----------------------------------------------------------------------------

Base.@kwdef mutable struct NEMODefault <: AbstractClosure
    Uv::Float64 = 0.1          # rn_Uv (m/s)
    Ue::Float64 = 0.02         # rn_Ue (m/s)
    Le::Union{Nothing, Float64} = nothing   # rn_Le (m); nothing → Δx
    gm::Bool = true
    ahm::Float64 = NaN         # filled by init_closure!
    aei::Float64 = NaN
end

default_stepper(::NEMODefault) = "RK4"

function init_closure!(c::NEMODefault, m)
    c.ahm = c.Uv * m.Δx^3 / 12
    c.aei = c.gm ? 0.5 * c.Ue * (isnothing(c.Le) ? m.Δx : c.Le) : 0.0
    return nothing
end

function pv_tendency!(Fqh, c::NEMODefault, f, m)
    # Constant coefficients → exact spectral forms.
    # bilaplacian on vorticity: −ahm ∇⁴ζ
    @. Fqh = -c.ahm * m.grid.Krsq^2 * f.ζh
    if c.gm
        Δψh = f.ψh[:, :, 1] .- f.ψh[:, :, 2]
        D = @. -c.aei * m.grid.Krsq * Δψh     # ∇·(κ∇Δψ) for constant κ
        @views @. Fqh[:, :, 1] += -m.F[1] * D
        @views @. Fqh[:, :, 2] +=  m.F[2] * D
    end
    Fqh .*= m.mask
    return Fqh
end

function param_pv_flux(c::NEMODefault, f, m)
    c.gm || return (0.0, 0.0)
    _, gy = interface_gradient(f.ψh[:, :, 1] .- f.ψh[:, :, 2], m)
    return gm_param_flux(c.aei, gy, m)
end

closure_diagnostics(c::NEMODefault, m) = (ahm = c.ahm, kappa_gm = c.aei)
