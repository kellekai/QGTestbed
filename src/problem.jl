# ----------------------------------------------------------------------------
# Problem construction
# ----------------------------------------------------------------------------

"""
    make_closure(name::Symbol; kwargs...)

`:none`, `:nemo`, `:geometric`, `:geometric_nosplit`, `:meke`, `:reynolds`.
"""
function make_closure(name::Symbol; kwargs...)
    name === :none              && return NoClosure()
    name === :nemo              && return NEMODefault(; kwargs...)
    name === :geometric         && return GEOMETRIC(; splitting = true, kwargs...)
    name === :geometric_nosplit && return GEOMETRIC(; splitting = false, kwargs...)
    name === :meke              && return MEKEGMBackscatter(; kwargs...)
    name === :reynolds          && return ReynoldsBackscatter(; kwargs...)
    error("unknown closure $name")
end

"Time step from an advective CFL estimate, capped at 1 h."
function default_dt(pp, nx; cfl = 0.25)
    Δx = pp.L / nx
    Uscale = max(0.2, 4 * pp.ΔU)
    dt = min(3600.0, cfl * Δx / Uscale)
    return 3600.0 / ceil(3600.0 / dt)       # divides an hour evenly
end

"""
    build_problem(pp, nx; closure = NoClosure(), dev = CPU(), dt = default_dt(pp, nx),
                  stepper = default_stepper(closure), T = Float64)

Two-layer MultiLayerQG problem on a doubly periodic β-plane with imposed
zonal shear `pp.U`, linear bottom drag `pp.μ`, and the closure attached
through the PV forcing hook `calcFq`. Returns `(prob, ctx, closure)`.
"""
function build_problem(pp, nx; closure::AbstractClosure = NoClosure(), dev = CPU(),
                       dt = default_dt(pp, nx), stepper = default_stepper(closure),
                       T = Float64)
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
                                μ = pp.μ, dt, stepper, T, forcing_kw...)

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
