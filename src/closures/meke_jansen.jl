# ----------------------------------------------------------------------------
# MEKE with GM + backscatter (Jansen, Adcroft, Khani & Kong 2019, JAMES),
# adapted to 2 layers.
#
# e : vertically averaged subgrid eddy kinetic energy (m² s⁻²)
#
#   K_GM   = c √(2e) L_mix R(Δ k_d)            GM diffusivity
#   ν₂     = −c √(2e) L_mix                     negative Laplacian viscosity (backscatter),
#                                               NOT scaled by R (as in the authors' MOM6 code)
#   ν₄     = max(c_smag |D|, 1/τ_vis) Δ⁴        biharmonic Smagorinsky (MOM6: max of the two)
#   L_mix  = min(Δ, (2e)^{1/4} β^{-1/2})        (flat bottom: β* = β)
#   R      = [1 + π/(Δ k_d)]⁻¹,  k_d = 1/Ld
#
#   ∂e/∂t = ė_GM + ė_ν4 − ė_ν2 − ė_diss − ∇·(u_bt e − K_e ∇e)
#     ė_GM   = K_GM f₀² |∇Δψ|² / (g′ H)          (PE release, vertical mean)
#     ė_ν4   = Σ (H_j/H) ν₄ (|∇²u_j|²)           (all biharmonic dissipation)
#     ė_ν2   = Σ (H_j/H) (−ν₂) |∇u_j|²           (energy given back by backscatter)
#     ė_diss = (c_d / H) √(2e + |u_bot|² + U_bg²) e
#     K_e    = c √(2e) L_mix
#
# Settings follow the authors' published MOM6 parameter files (github.com/MFJansen/
# MEKE_GM_BS_parameters, e.g. MOM_parameter_doc.all_p25_meke_GM_BS) and their MOM6 code
# (github.com/MFJansen/MOM6), which produced the paper's results:
#   MEKE_KHCOEFF = 0.15, MEKE_VISCOSITY_COEFF = −0.15  → c = 0.15 for K_GM and ν₂
#   KH_RES_SCALE_COEF = π, KH_RES_FN_POWER = 1, RESOLN_SCALED_KHTH = True → R on K_GM
#   MOM_hor_visc adds MEKE%Ku to the viscosity without the resolution function → ν₂ unscaled
#   SMAG_BI_CONST = 0.06, AH_TIME_SCALE = 45 d, Ah = MAX(background, Smagorinsky)
#   CDRAG = 0.003, DRAG_BG_VEL = 0.1 m/s, MEKE_CD_SCALE = 1, MEKE_USCALE = 0:
#     drag rate = CDRAG/H · √(2·MEKE + u_bot² + U_bg²)   (MOM_MEKE.F90; 1/H via I_mass)
#   MEKE_MIN_LSCALE = True, MEKE_ALPHA_GRID = MEKE_ALPHA_RHINES = 1 → L_mix = min(Δ, L_Rhines)
#   MEKE_CT = MEKE_CB = 0 (no barotropic/bottom factors), MEKE_KHMEKE_FAC = 1 (K_e = c√(2e)L_mix)
# Differences to the paper's text, which states c_d = 0.03, ν₂ = −K_GM (with R), √(e + …) and
# ν₄ = (c_smag|D| + 1/τ)Δ⁴: the parameter files and code are what was run, so they are the
# defaults here. The paper-text variant: c_d = 0.03, bs_resolution_scaled = true,
# drag_eke_factor = 1, nu4_max = false.
# Backscatter is applied with the same ν₂ in both
# layers (the paper does not prescribe a vertical structure; Bachman 2019
# used barotropic-only, available via `bs_barotropic_only = true`).
# ----------------------------------------------------------------------------

Base.@kwdef mutable struct MEKEGMBackscatter <: AbstractClosure
    c::Float64 = 0.15
    c_smag::Float64 = 0.06
    τ_vis::Float64 = 45 * 86400.0
    c_d::Float64 = 0.003               # CDRAG of the authors' runs (paper text: 0.03)
    U_bg::Float64 = 0.1                # DRAG_BG_VEL
    e_init::Float64 = 1e-4
    resolution_function::Bool = true
    bs_resolution_scaled::Bool = false # true: ν₂ = −K_GM incl. R (paper text); false: MOM6 code
    drag_eke_factor::Float64 = 2.0     # √(f·e + …): 2 in MOM6 (2·MEKE_CD_SCALE²·MEKE), 1 in the paper text
    nu4_max::Bool = true               # ν₄ = max(c_smag|D|, 1/τ)Δ⁴ (MOM6); false: sum (paper text)
    bs_barotropic_only::Bool = false
    e::Any = nothing
    K::Any = nothing
    budget::Any = nothing      # domain means of the MEKE budget terms at the last step (diagnostics)
end

# Names of the MEKE budget diagnostics (domain means, written to the time series):
#   meke_src_GM        ė_GM   PE released by GM                       (m² s⁻³, source)
#   meke_src_biharm    ė_ν4   energy removed by biharmonic Smagorinsky (m² s⁻³, source)
#   meke_sink_bs       ė_ν2   energy given back by backscatter         (m² s⁻³, sink)
#   meke_sink_drag     ė_diss bottom-drag dissipation of MEKE          (m² s⁻³, sink)
#   meke_net           sum of the four (≈ d⟨e⟩/dt; domain-mean transport is zero)
#   meke_drag_rate     c_d/H √(2e + |u_bot|² + U_bg²)                  (s⁻¹)
#   meke_drag_rate_Ubg c_d/H · U_bg, the part set by U_bg alone        (s⁻¹)
#   meke_Lmix          mixing length L_mix                              (m)
const MEKE_BUDGET_KEYS = (:meke_src_GM, :meke_src_biharm, :meke_sink_bs, :meke_sink_drag, :meke_net,
                          :meke_drag_rate, :meke_drag_rate_Ubg, :meke_Lmix)

default_stepper(::MEKEGMBackscatter) = "RK4"
has_prognostic_state(::MEKEGMBackscatter) = true

function init_closure!(c::MEKEGMBackscatter, m)
    c.e = fill!(similar(m.real2), c.e_init)
    c.K = fill!(similar(m.real2), 0)
    c.budget = NamedTuple{MEKE_BUDGET_KEYS}(ntuple(_ -> NaN, length(MEKE_BUDGET_KEYS)))
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

# Velocity derivatives are taken directly from the spectral velocities of `layerfields` (one inverse
# FFT each, no physical→spectral round trip) and shared by the Smagorinsky viscosity, the biharmonic
# term and the backscatter. Same result as the round-trip form up to roundoff, about half the FFTs.

"Velocity gradients (physical, (nx, ny, 2)) from the spectral layer velocities."
function velocity_gradients(f, m)
    kr, l = m.grid.kr, m.grid.l
    return (; ux = tophys(@.(im * kr * f.uh), m), uy = tophys(@.(im * l * f.uh), m),
              vx = tophys(@.(im * kr * f.vh), m), vy = tophys(@.(im * l * f.vh), m))
end

function smagorinsky_ν4(c::MEKEGMBackscatter, g, m)
    D = @. sqrt((g.ux - g.vy)^2 + (g.uy + g.vx)^2)
    c.nu4_max && return @. max(c.c_smag * D, 1 / c.τ_vis) * m.Δx^4
    return @. (c.c_smag * D + 1 / c.τ_vis) * m.Δx^4
end

"Backscatter viscosity magnitude |ν₂| = c √(2e) L_mix (× R only if bs_resolution_scaled)."
function meke_Kbs(c::MEKEGMBackscatter, m)
    c.bs_resolution_scaled && return meke_K(c, m)
    Lmix = meke_Lmix(c, m)
    return @. c.c * sqrt(2 * c.e) * Lmix
end

"Gradients the backscatter acts on: per layer, or of the barotropic velocity in both layers."
function backscatter_gradients(c::MEKEGMBackscatter, g, m)
    c.bs_barotropic_only || return g
    bt(A) = (a = layer_average(A, m); cat(a, a; dims = 3))
    return (; ux = bt(g.ux), uy = bt(g.uy), vx = bt(g.vx), vy = bt(g.vy))
end

"""
MEKE closure terms on the state `f`: K_GM, interface gradient, biharmonic Smagorinsky and
negative-Laplacian backscatter with their dissipation densities, and (if `tendencies`) their
spectral vorticity tendencies. step_closure! needs only the dissipations (tendencies = false).
"""
function meke_terms(c::MEKEGMBackscatter, f, m; tendencies = true)
    K = meke_K(c, m)
    Δψh = f.ψh[:, :, 1] .- f.ψh[:, :, 2]
    gx, gy = interface_gradient(Δψh, m)
    g = velocity_gradients(f, m)
    ν4 = smagorinsky_ν4(c, g, m)
    lu = tophys(@.(-m.grid.Krsq * f.uh), m)                         # ∇²u, ∇²v
    lv = tophys(@.(-m.grid.Krsq * f.vh), m)
    diss4 = @. ν4 * (lu^2 + lv^2)
    b = backscatter_gradients(c, g, m)
    ν2 = -meke_Kbs(c, m)                                            # ν₂ = −c√(2e)L_mix
    diss2 = @. ν2 * (b.ux^2 + b.uy^2 + b.vx^2 + b.vy^2)
    tendencies || return (; K, gx, gy, diss4, diss2)
    kr, l = m.grid.kr, m.grid.l
    tu4 = @. m.grid.Krsq * $toh(ν4 .* lu, m)                        # spectral −∇²(ν4∇²u)
    tv4 = @. m.grid.Krsq * $toh(ν4 .* lv, m)
    Fζ4h = @. im * kr * tv4 - im * l * tu4
    tu2 = divh(ν2 .* b.ux, ν2 .* b.uy, m)                           # spectral ∇·(ν₂∇u)
    tv2 = divh(ν2 .* b.vx, ν2 .* b.vy, m)
    Fζ2h = @. im * kr * tv2 - im * l * tu2
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
    t = meke_terms(c, f, m; tendencies = false)
    c.K .= t.K
    ė_GM  = gm_energy_release(t.K, t.gx, t.gy, m) ./ m.Htot
    ė_ν4  = layer_average(t.diss4, m)
    ė_ν2  = -layer_average(t.diss2, m)          # diss2 < 0 for ν₂ < 0 → ė_ν2 > 0
    ubot2 = @views @. (f.u[:, :, 2] + m.U[2])^2 + f.v[:, :, 2]^2
    drag_rate = @. c.c_d / m.Htot * sqrt(c.drag_eke_factor * c.e + ubot2 + c.U_bg^2)
    ė_diss = @. drag_rate * c.e
    source = @. ė_GM + ė_ν4 - ė_ν2 - ė_diss
    c.budget = NamedTuple{MEKE_BUDGET_KEYS}((mean(ė_GM), mean(ė_ν4), mean(ė_ν2), mean(ė_diss), mean(source),
                                             mean(drag_rate), c.c_d / m.Htot * c.U_bg, mean(meke_Lmix(c, m))))
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

closure_state(c::MEKEGMBackscatter) = Array(c.e)
set_closure_state!(c::MEKEGMBackscatter, s) = (copyto!(c.e, s); nothing)

closure_diagnostics(c::MEKEGMBackscatter, m) =
    merge((MEKE = mean(c.e), kappa_gm = mean(c.K), R = meke_R(c, m)), c.budget)
