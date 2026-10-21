# ----------------------------------------------------------------------------
# GEOMETRIC (Mak et al. 2023, JAMES, arXiv:2306.08988), adapted to 2 layers
#
#   κ_gm = α Ê / ∫(M²/N) dz
#   ∂Ê/∂t + ∇·[(ū_bt − |c| x̂) Ê] = ∫ κ_gm M⁴/N² dz − λ(Ê − Ê₀) + η_E ∇²Ê
#
# Ê is the depth-integrated parametrized eddy energy (m³ s⁻²).
#
# 2-layer translation (N Δz = √(g′H), M²/N = f₀|Δu|/√(g′H), Δu = u1 − u2):
#   ∫(M²/N) dz      = f₀ |Δu| √(H/g′)
#   κ_gm            = α Ê √(g′/H) / (f₀ |Δu|)
#   ∫κ M⁴/N² dz     = κ f₀² |Δu|² / g′   (= exact PE loss of the resolved flow)
#
# Scale awareness ("splitting", Mak et al. 2023, Sec. 3):
#   the stratification is split into a large-scale part with
#   (1 − L²∇²)^M Θ_L = Θ, L = 100 km, M = 2. GM uses the large-scale slope
#   only, in both κ_gm and the energy source, so resolved eddies are left alone.
#   In QG the "stratification" is Δψ; the filter is exact in spectral space:
#   Δψ̂_L = Δψ̂ / (1 + L²k²)^M. The background shear is large scale by definition.
#   `splitting = false` gives the original GEOMETRIC.
#
# Defaults follow the authors' NEMO 4.0.5 code and namelists for the eddy-permitting runs
# (MY_SRC/ldfeke.F90, ldftra_split.F90; EXP00/namelist_cfg_R010/R025/R050):
#   rn_geom = 0.06 (α), rn_ekedis = 80 d (λ⁻¹), rn_eke_lap = 500 m²/s (η_E),
#   rn_eke_min = 5 m³/s² (Ê₀; the paper's table gives 4), nn_eke_opt = 1 (PE→EKE, dissipation,
#   advection by the depth-mean flow and the long Rossby speed; no KE→EKE source),
#   rn_aeiv_min = 10 m²/s, rn_aeiv_max = 3·10⁴ m²/s (bounds on κ_gm),
#   splitting: (1 − L²∇²) solved twice = (1 + L²k²)⁻² on T and S, rn_L = 100 km (R025, R050;
#   200 km at R010), slopes and N² for κ_gm, the source and the GM velocities from the large scale,
#   momentum: bilaplacian, nn_ahm_ijk_t = 0, ahm = U_v Δx³/12 with rn_Uv = 0.5 m/s, rn_Lv = Δx.
# Not reproduced: the taper of κ_gm where 0.5∫N dz/f < rn_ross_min = 7 km (that radius is
# ≥ 19 km at every testbed latitude, so the taper is 1).
# ----------------------------------------------------------------------------

Base.@kwdef mutable struct GEOMETRIC <: AbstractClosure
    α::Float64 = 0.06
    λ::Float64 = 1 / (80 * 86400.0)
    ηE::Float64 = 500.0
    E0::Float64 = 5.0          # rn_eke_min of the authors' runs (paper table: 4)
    splitting::Bool = true
    L_split::Float64 = 100e3
    M_split::Int = 2
    Uv::Float64 = 0.5          # rn_Uv of the authors' eddy-permitting runs (bilaplacian, Lv = Δx)
    κ_max::Float64 = 3.0e4     # rn_aeiv_max
    κ_min::Float64 = 10.0      # rn_aeiv_min
    shear_min::Float64 = 1e-4  # m/s, avoids division by zero in κ
    westward_advection::Bool = true
    ahm::Float64 = NaN
    E::Any = nothing           # Ê, (nx, ny), on the model device
    κ::Any = nothing           # last κ_gm, for diagnostics
end

default_stepper(::GEOMETRIC) = "RK4"
has_prognostic_state(::GEOMETRIC) = true

function init_closure!(c::GEOMETRIC, m)
    c.ahm = c.Uv * m.Δx^3 / 12
    c.E = fill!(similar(m.real2), c.E0)
    c.κ = fill!(similar(m.real2), 0)
    return nothing
end

function geometric_split(Δψh, c::GEOMETRIC, m)
    c.splitting || return Δψh
    return @. Δψh / (1 + c.L_split^2 * m.grid.Krsq)^c.M_split
end

function geometric_kappa(gx, gy, c::GEOMETRIC, m)
    shear = @. max(sqrt(gx^2 + gy^2), c.shear_min)
    κ = @. c.α * c.E * sqrt(m.g′ / m.Htot) / (m.f₀ * shear)
    return clamp.(κ, c.κ_min, c.κ_max)
end

function pv_tendency!(Fqh, c::GEOMETRIC, f, m)
    @. Fqh = -c.ahm * m.grid.Krsq^2 * f.ζh                 # bilaplacian (constant)
    ΔψLh = geometric_split(f.ψh[:, :, 1] .- f.ψh[:, :, 2], c, m)
    gx, gy = interface_gradient(ΔψLh, m)
    κ = geometric_kappa(gx, gy, c, m)
    gm_pv_tendency!(Fqh, κ, gx, gy, m)
    Fqh .*= m.mask
    return Fqh
end

function step_closure!(c::GEOMETRIC, f, m, dt)
    ΔψLh = geometric_split(f.ψh[:, :, 1] .- f.ψh[:, :, 2], c, m)
    gx, gy = interface_gradient(ΔψLh, m)
    κ = geometric_kappa(gx, gy, c, m)
    c.κ .= κ
    source = gm_energy_release(κ, gx, gy, m) .- c.λ .* (c.E .- c.E0)
    ubt, vbt = barotropic_velocity(f, m)
    if c.westward_advection
        ubt = ubt .- m.pp.c_rossby
    end
    tend_h = advect_diffuse(c.E, ubt, vbt, c.ηE, m)
    euler_step_positive!(c.E, tend_h, source, dt, m)
    return nothing
end

function param_pv_flux(c::GEOMETRIC, f, m)
    ΔψLh = geometric_split(f.ψh[:, :, 1] .- f.ψh[:, :, 2], c, m)
    gx, gy = interface_gradient(ΔψLh, m)
    return gm_param_flux(geometric_kappa(gx, gy, c, m), gy, m)
end

closure_state(c::GEOMETRIC) = Array(c.E)
set_closure_state!(c::GEOMETRIC, s) = (copyto!(c.E, s); nothing)

closure_diagnostics(c::GEOMETRIC, m) =
    (E_param = mean(c.E), kappa_gm = mean(c.κ), kappa_gm_max = maximum(c.κ))
