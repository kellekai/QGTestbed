# ----------------------------------------------------------------------------
# Closure interface
#
#   pv_tendency!(Fqh, closure, f, m)  : write the PV tendency (spectral,
#                                        (nkr, nl, 2)) into Fqh. Called inside
#                                        every RK stage; must NOT change state.
#   step_closure!(closure, f, m, dt)  : advance prognostic closure state
#                                        (eddy energy) once per model step.
#   init_closure!(closure, m)         : allocate state on the model device.
#   default_stepper(closure)          : time stepper / base dissipation.
#   closure_diagnostics(closure, m)   : NamedTuple of scalars to log.
#
# `f` is the output of `layerfields`, `m` a `ModelContext`.
# ----------------------------------------------------------------------------

abstract type AbstractClosure end

"Reference: no closure, only the scale-selective spectral filter of FilteredRK4."
struct NoClosure <: AbstractClosure end

pv_tendency!(Fqh, ::AbstractClosure, f, m) = (Fqh .= 0; Fqh)
step_closure!(::AbstractClosure, f, m, dt) = nothing
init_closure!(::AbstractClosure, m) = nothing
closure_diagnostics(::AbstractClosure, m) = NamedTuple()

# Closures that bring their own dissipation run with plain RK4; the others use
# the exponential spectral filter (as pyqg does, and as in Perezhogin & Glazunov 2023).
default_stepper(::AbstractClosure) = "FilteredRK4"

"Closures with prognostic eddy energy need fields every step."
has_prognostic_state(::AbstractClosure) = false

"""
    closure_state(closure) / set_closure_state!(closure, state)

The prognostic eddy-energy field (CPU array), stored by `run_simulation` next to every PV
snapshot so the closure tendency can be recomputed offline (e.g. its spectral energy
transfer). `nothing` for closures without prognostic state.
"""
closure_state(::AbstractClosure) = nothing
set_closure_state!(::AbstractClosure, state) = nothing

"""
    param_pv_flux(closure, f, m) -> (flux_layer1, flux_layer2)

Domain-mean meridional PV flux carried by the closure (zero for closures
without one). With the imposed shear held fixed, the GM part of this flux acts
on the mean state that the model never updates, so it is invisible in the
resolved fields; compare resolved ⟨v q⟩ + this flux with the truth's ⟨v q⟩.
Momentum (bi)harmonic viscosity terms are not included.
"""
param_pv_flux(::AbstractClosure, f, m) = (0.0, 0.0)

"""
GM PV flux per layer, from ∂q₁ = −∇·(F₁κ∇Δψ), ∂q₂ = −∇·(−F₂κ∇Δψ):
meridional means ⟨F₁κ ∂yΔψ⟩ and ⟨−F₂κ ∂yΔψ⟩ with ∂yΔψ including −ΔU.
"""
gm_param_flux(κ, gy, m) = (mean(m.F[1] .* κ .* gy), mean(-m.F[2] .* κ .* gy))

"Scalar settings of a closure (for the output file)."
closure_settings(c::AbstractClosure) =
    Dict(String(k) => getfield(c, k) for k in fieldnames(typeof(c))
         if getfield(c, k) isa Union{Number, Bool, Symbol, Nothing})

# ------------------------------------------------------------------------------
# Building blocks shared by several closures
# ------------------------------------------------------------------------------

"""
    interface_gradient(Δψh, m)

Gradient of Δψ_total = ψ1 − ψ2 including the imposed shear
(Δψ_bg = −ΔU y), so |∇Δψ_total| = |u1 − u2| is the local vertical shear.
Returns physical (gx, gy).
"""
function interface_gradient(Δψh, m)
    gx = tophys(@.(im * m.grid.kr * Δψh), m)
    gy = tophys(@.(im * m.grid.l  * Δψh), m) .- m.ΔU
    return gx, gy
end

"""
    gm_pv_tendency!(Fqh, κ, gx, gy, m)

QG analogue of Gent–McWilliams: down-gradient diffusion of interface height,
i.e. of the stretching part of PV. With Δψ = ψ1 − ψ2,

    ∂q1/∂t = −F1 ∇·(κ∇Δψ),   ∂q2/∂t = +F2 ∇·(κ∇Δψ).

The background gradient is included in (gx, gy); its flux divergence only
matters where κ varies. The resolved-energy loss is exactly
∫ κ f₀² |∇Δψ|² / g′ dA (depth integrated), which is what the eddy-energy
budgets receive as source.
"""
function gm_pv_tendency!(Fqh, κ, gx, gy, m)
    D = divh(κ .* gx, κ .* gy, m)              # 2D spectral ∇·(κ∇Δψ)
    @views @. Fqh[:, :, 1] += -m.F[1] * D
    @views @. Fqh[:, :, 2] +=  m.F[2] * D
    return Fqh
end

"Depth-integrated PE release by GM: κ f₀² |∇Δψ|² / g′ (m³ s⁻³)."
gm_energy_release(κ, gx, gy, m) = @. κ * m.f₀^2 * (gx^2 + gy^2) / m.g′

"Add a spectral vorticity tendency (same for both layers if 2D) to Fqh."
function add_vorticity_tendency!(Fqh, Fζh)
    Fqh .+= Fζh
    return Fqh
end

"""
    advect_diffuse(E, ux, uy, ηE, m)

Spectral tendency −∇·(u E) + ∇·(ηE ∇E) for a 2D field E with 2D velocity
(ux, uy) and scalar or 2D diffusivity ηE.
"""
function advect_diffuse(E, ux, uy, ηE, m)
    Ex, Ey = ∂x(E, m), ∂y(E, m)
    return divh(@.(ηE * Ex - ux * E), @.(ηE * Ey - uy * E), m)
end

"Barotropic (thickness-weighted) velocity including the imposed flow."
function barotropic_velocity(f, m)
    ubt = layer_average(f.u, m) .+ (m.H[1] * m.U[1] + m.H[2] * m.U[2]) / m.Htot
    vbt = layer_average(f.v, m)
    return ubt, vbt
end

"Forward-Euler step of a positive 2D field, applied in physical space."
function euler_step_positive!(E, tendency_h, source_phys, dt, m; Emin = 0)
    E .+= dt .* (tophys(m.mask .* tendency_h, m) .+ source_phys)
    E .= max.(E, Emin)
    return E
end
