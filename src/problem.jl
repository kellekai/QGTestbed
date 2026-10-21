# ----------------------------------------------------------------------------
# Problem construction
# ----------------------------------------------------------------------------

"""
    make_closure(name::Symbol; kwargs...)

`:none`, `:nemo`, `:nemo_eorca025`, `:geometric`, `:geometric_nosplit`, `:meke`, `:reynolds`.
`:nemo_eorca025` is NEMO as in the project's eORCA025 runs (namelist a3ue): no GM
(ln_ldfeiv = .false.), bilaplacian momentum viscosity with nn_ahm_ijk_t = 20, rn_Uv = 0.0838 m/s.
`:geometric_eorca025` is GEOMETRIC (Mak et al. 2023, with splitting) with that eORCA025 viscosity
instead of the authors' rn_Uv = 0.5.
"""
function make_closure(name::Symbol; kwargs...)
    name === :none              && return NoClosure()
    name === :nemo              && return NEMODefault(; kwargs...)
    name === :nemo_eorca025     && return NEMODefault(; gm = false, Uv = 0.0838, kwargs...)
    name === :geometric         && return GEOMETRIC(; splitting = true, kwargs...)
    name === :geometric_nosplit && return GEOMETRIC(; splitting = false, kwargs...)
    name === :geometric_eorca025 && return GEOMETRIC(; splitting = true, Uv = 0.0838, kwargs...)
    name === :meke              && return MEKEGMBackscatter(; kwargs...)
    name === :reynolds          && return ReynoldsBackscatter(; kwargs...)
    error("unknown closure $name")
end

"""
    default_aliased_fraction()

Dealiasing of the nonlinear term (GeophysicalFlows `aliased_fraction`). Default 0: no 2/3-rule
dealiasing, only the exponential filter, as in pyqg. This setting reproduces pyqg and
Perezhogin & Glazunov (2023, Fig. 10) at 64² (scripts/perezhogin_validation.sh, pyqg-validation/);
the 2/3 rule (1/3) made coarse runs too dissipative. Override with the environment variable
QG_ALIASED_FRACTION (e.g. 0.3333333 to reproduce runs made before this change).
"""
default_aliased_fraction() = parse(Float64, get(ENV, "QG_ALIASED_FRACTION", "0"))

"""
    filter_all()

`true` when the environment variable QG_FILTER_ALL is true/1/yes: every closure then runs with
the exponential filter (FilteredRK4), as all runs in pyqg do, also those that bring their own
dissipation (nemo, geometric*, meke). Default false: those closures run with plain RK4.
"""
filter_all() = lowercase(get(ENV, "QG_FILTER_ALL", "false")) in ("1", "true", "yes")

"Time stepper of a closure, honouring QG_FILTER_ALL."
choose_stepper(closure) = filter_all() ? "FilteredRK4" : default_stepper(closure)

"Time step from an advective CFL estimate, capped at 1 h."
function default_dt(pp, nx; cfl = 0.25)
    Δx = pp.L / nx
    Uscale = max(0.2, 4 * pp.ΔU)
    dt = min(3600.0, cfl * Δx / Uscale)
    return 3600.0 / ceil(3600.0 / dt)       # divides an hour evenly
end

"""
    build_problem(pp, nx; closure = NoClosure(), dev = CPU(), dt = default_dt(pp, nx),
                  stepper = choose_stepper(closure), T = Float64,
                  aliased_fraction = default_aliased_fraction())

Two-layer MultiLayerQG problem on a doubly periodic β-plane with imposed
zonal shear `pp.U`, linear bottom drag `pp.μ`, and the closure attached
through the PV forcing hook `calcFq`. Returns `(prob, ctx, closure)`.
`aliased_fraction` is GeophysicalFlows' dealiasing: 0 (default, see `default_aliased_fraction`)
= none, as in pyqg; 1/3 = 2/3 rule. `stepper` follows `choose_stepper` (QG_FILTER_ALL).
"""
function build_problem(pp, nx; closure::AbstractClosure = NoClosure(), dev = CPU(),
                       dt = default_dt(pp, nx), stepper = choose_stepper(closure),
                       T = Float64, aliased_fraction = default_aliased_fraction())
    ctx_ref = Ref{Any}(nothing)

    # Without a closure, leave calcFq at the package default (no forcing arrays).
    forcing_kw = if closure isa NoClosure
        NamedTuple()
    else
        calcFq! = function (Fqh, sol, t, clock, vars, params, grid)
            m = ctx_ref[]
            pv_tendency!(Fqh, closure, layerfields(sol, m), m)
            return nothing
        end
        (; calcFq = calcFq!)
    end

    prob = MultiLayerQG.Problem(2, dev;
                                nx, Lx = pp.L, ny = nx, Ly = pp.L,
                                f₀ = pp.f₀, β = pp.β,
                                U = pp.U, H = [pp.H1, pp.H2],
                                b = [0.0, -pp.g′],      # g′ = b₁ − b₂
                                μ = pp.μ, dt, stepper, T, aliased_fraction, forcing_kw...)

    ctx = ModelContext(prob, pp)
    ctx_ref[] = ctx
    init_closure!(closure, ctx)
    return prob, ctx, closure
end

"Random small PV perturbation (as in pyqg), reproducible via `seed`."
function set_random_initial_condition!(prob, pp; amplitude = 1e-7, seed = 1)
    rng = MersenneTwister(seed)
    nx, ny = prob.grid.nx, prob.grid.ny
    q0 = amplitude .* randn(rng, eltype(prob.grid), nx, ny, 2)
    q0 .-= mean(q0; dims = (1, 2))
    MultiLayerQG.set_q!(prob, device_array(prob.grid.device)(q0))
    return prob
end
