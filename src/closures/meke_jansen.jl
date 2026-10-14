# ----------------------------------------------------------------------------
# MEKE with GM + backscatter (Jansen, Adcroft, Khani & Kong 2019, JAMES),
# adapted to 2 layers.
#
# e : vertically averaged subgrid eddy kinetic energy (m² s⁻²)
#
#   K_GM   = c √(2e) L_mix R(Δ k_d)            GM diffusivity
#   ν₂     = −K_GM                              negative Laplacian viscosity (backscatter)
#   ν₄     = (c_smag |D| + 1/τ_vis) Δ⁴          biharmonic Smagorinsky
#   L_mix  = min(Δ, (2e)^{1/4} β^{-1/2})        (flat bottom: β* = β)
#   R      = [1 + π/(Δ k_d)]⁻¹,  k_d = 1/Ld
#
#   ∂e/∂t = ė_GM + ė_ν4 − ė_ν2 − ė_diss − ∇·(u_bt e − K_e ∇e)
#     ė_GM   = K_GM f₀² |∇Δψ|² / (g′ H)          (PE release, vertical mean)
#     ė_ν4   = Σ (H_j/H) ν₄ (|∇²u_j|²)           (all biharmonic dissipation)
#     ė_ν2   = Σ (H_j/H) (−ν₂) |∇u_j|²           (energy given back by backscatter)
#     ė_diss = (c_d / H) √(e + |u_bot|² + U_bg²) e
#     K_e    = c √(2e) L_mix
#
# Paper values: c = 0.15, c_smag = 0.06, τ_vis = 45 d, U_bg = 0.1 m/s,
# c_d = 0.03 (stated in the paper's text after its Eq. 9). The paper writes
# ė_diss without a depth factor; e is a vertical mean, so 1/H is needed for
# consistent units (as in the MOM6 MEKE implementation) and is applied here.
# Backscatter is applied with the same ν₂ in both
# layers (the paper does not prescribe a vertical structure; Bachman 2019
# used barotropic-only, available via `bs_barotropic_only = true`).
# ----------------------------------------------------------------------------

Base.@kwdef mutable struct MEKEGMBackscatter <: AbstractClosure
    c::Float64 = 0.15
    c_smag::Float64 = 0.06
    τ_vis::Float64 = 45 * 86400.0
    c_d::Float64 = 0.03
    U_bg::Float64 = 0.1
    e_init::Float64 = 1e-4
    resolution_function::Bool = true
    bs_barotropic_only::Bool = false
    e::Any = nothing
    K::Any = nothing
end

default_stepper(::MEKEGMBackscatter) = "RK4"
has_prognostic_state(::MEKEGMBackscatter) = true

function init_closure!(c::MEKEGMBackscatter, m)
    c.e = fill!(similar(m.real2), c.e_init)
    c.K = fill!(similar(m.real2), 0)
    return nothing
end

meke_R(c::MEKEGMBackscatter, m) = c.resolution_function ? 1 / (1 + π * m.Ld / m.Δx) : 1.0

function meke_Lmix(c::MEKEGMBackscatter, m)
    return @. min(m.Δx, (2 * c.e)^(1 / 4) / sqrt(m.β))
end

function meke_K(c::MEKEGMBackscatter, m)
    Lmix = meke_Lmix(c, m)
    R = meke_R(c, m)
    return @. c.c * sqrt(2 * c.e) * Lmix * R
end

function smagorinsky_ν4(c::MEKEGMBackscatter, f, m)
    ux, uy, vx, vy = ∂x(f.u, m), ∂y(f.u, m), ∂x(f.v, m), ∂y(f.v, m)
    D = @. sqrt((ux - vy)^2 + (uy + vx)^2)
    return @. (c.c_smag * D + 1 / c.τ_vis) * m.Δx^4
end

"Backscatter velocities: per layer, or the barotropic part only."
function backscatter_velocities(c::MEKEGMBackscatter, f, m)
    c.bs_barotropic_only || return f.u, f.v
    ubt, vbt = layer_average(f.u, m), layer_average(f.v, m)
    return cat(ubt, ubt; dims = 3), cat(vbt, vbt; dims = 3)
end

function meke_terms(c::MEKEGMBackscatter, f, m)
    K = meke_K(c, m)
    Δψh = f.ψh[:, :, 1] .- f.ψh[:, :, 2]
    gx, gy = interface_gradient(Δψh, m)
    ν4 = smagorinsky_ν4(c, f, m)
    Fζ4h, diss4 = curl_of_momentum_biharmonic(ν4, f.u, f.v, m)
    ub, vb = backscatter_velocities(c, f, m)
    Fζ2h, diss2 = curl_of_momentum_laplacian(-K, ub, vb, m)   # ν₂ = −K_GM
    return (; K, gx, gy, Fζ4h, diss4, Fζ2h, diss2)
end

function pv_tendency!(Fqh, c::MEKEGMBackscatter, f, m)
    t = meke_terms(c, f, m)
    Fqh .= t.Fζ4h .+ t.Fζ2h
    gm_pv_tendency!(Fqh, t.K, t.gx, t.gy, m)
    Fqh .*= m.mask
    return Fqh
end

function step_closure!(c::MEKEGMBackscatter, f, m, dt)
    t = meke_terms(c, f, m)
    c.K .= t.K
    ė_GM  = gm_energy_release(t.K, t.gx, t.gy, m) ./ m.Htot
    ė_ν4  = layer_average(t.diss4, m)
    ė_ν2  = -layer_average(t.diss2, m)          # diss2 < 0 for ν₂ < 0 → ė_ν2 > 0
    ubot2 = @views @. (f.u[:, :, 2] + m.U[2])^2 + f.v[:, :, 2]^2
    ė_diss = @. c.c_d / m.Htot * sqrt(c.e + ubot2 + c.U_bg^2) * c.e
    source = @. ė_GM + ė_ν4 - ė_ν2 - ė_diss
    ubt, vbt = barotropic_velocity(f, m)
    Ke = @. c.c * sqrt(2 * c.e) * $meke_Lmix(c, m)
    tend_h = advect_diffuse(c.e, ubt, vbt, Ke, m)
    euler_step_positive!(c.e, tend_h, source, dt, m; Emin = 1e-10)
    return nothing
end

function param_pv_flux(c::MEKEGMBackscatter, f, m)
    Δψh = f.ψh[:, :, 1] .- f.ψh[:, :, 2]
    _, gy = interface_gradient(Δψh, m)
    return gm_param_flux(meke_K(c, m), gy, m)
end

closure_diagnostics(c::MEKEGMBackscatter, m) =
    (MEKE = mean(c.e), kappa_gm = mean(c.K), R = meke_R(c, m))
