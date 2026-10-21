# ----------------------------------------------------------------------------
# Driver for one sweep member (one latitude, one resolution, one closure)
# ----------------------------------------------------------------------------

const DAY = 86400.0
const YEAR = 365 * DAY

"""
    run_simulation(; lat = 49.0, nx = 64, closure = :none, years = 20.0,
                     spinup_years = 5.0, dev = CPU(), outdir = "output",
                     diag_every_days = 1.0, spectra_every_days = 10.0,
                     snapshot_every_days = Inf, vorticity_every_days = 30.0,
                     seed = 1, case_kwargs = (;), closure_kwargs = (;), dt = nothing)

Run one member and write `outdir/<lat>N_nx<nx>_<closure>_seed<seed>.jld2`
containing parameters, time series (energies, closure energy rate, closure
scalars, resolved and parametrized meridional PV flux per layer), time-mean
KE spectra and PV fluxes after spin-up, upper-layer relative vorticity
snapshots and lower-layer streamfunction snapshots (`psi_lower`) from t = 0
(Float32, every `vorticity_every_days`, including the spin-up), and (optional) full PV snapshots after spin-up for coarse-graining,
with the closure's eddy-energy field at the same times (`closure_snapshots`, geometric and meke only).
`closure = :truth` is a no-closure run at the resolution chosen by
`truth_resolution` when `nx = 0`. `tag` is appended to the file name and
`metadata` (e.g. tuning parameters) is stored in the file as given.
`aliased_fraction` (default 0 = no dealiasing, pyqg numerics; 1/3 = 2/3 rule, see
`default_aliased_fraction`), the stepper and whether the run stopped early (NaN) are stored too.
"""
function run_simulation(; lat = 49.0, nx = 64, closure = :none, years = 20.0,
                          spinup_years = 5.0, dev = CPU(), outdir = "output",
                          diag_every_days = 1.0, spectra_every_days = 10.0,
                          snapshot_every_days = Inf, vorticity_every_days = 30.0,
                          progress_every_days = 10.0,
                          seed = 1, case_kwargs = (;), closure_kwargs = (;), dt = nothing,
                          tag = "", metadata = Dict{String, Any}(),
                          aliased_fraction = default_aliased_fraction())
    pp = physical_parameters(CaseConfig(; lat, case_kwargs...))
    is_truth = closure === :truth
    if is_truth && nx == 0
        nx = truth_resolution(pp)
    end
    cl = make_closure(is_truth ? :none : closure; closure_kwargs...)
    dt = isnothing(dt) ? default_dt(pp, nx) : dt
    stepper = choose_stepper(cl)
    prob, m, cl = build_problem(pp, nx; closure = cl, dev, dt, aliased_fraction, stepper)
    set_random_initial_condition!(prob, pp; seed)

    name = @sprintf("%02.0fN_nx%d_%s_seed%d", lat, nx, String(closure), seed)
    isempty(tag) || (name *= "_" * tag)        # e.g. tuning runs: ..._seed1_C_R=7.0
    mkpath(outdir)
    file = joinpath(outdir, name * ".jld2")

    nsteps = round(Int, years * YEAR / dt)
    nspin = round(Int, spinup_years * YEAR / dt)
    nprog = max(1, round(Int, progress_every_days * DAY / dt))

    Δx = pp.L / nx
    log_line("="^78)
    log_line(@sprintf("Run %s", name))
    log_line(@sprintf("  latitude %.1f°N   f0 = %.3e s⁻¹   beta = %.3e m⁻¹s⁻¹", lat, pp.f₀, pp.β))
    log_line(@sprintf("  grid %d²   dx = %.2f km   Ld = %.2f km   dx/Ld = %.2f", nx, Δx / 1e3, pp.Ld / 1e3, Δx / pp.Ld))
    log_line(@sprintf("  shear dU = %.4f m/s   criticality xi = %.2f   dealiasing aliased_fraction = %.3f", pp.ΔU, pp.ξ, aliased_fraction))
    log_line(@sprintf("  closure %s   stepper %s   settings %s", String(closure), stepper,
                      join(["$k=$v" for (k, v) in closure_settings(cl)], ", ")))
    log_line(@sprintf("  dt = %.1f s   %d steps (%.1f yr, spin-up %.1f yr = %d steps)", dt, nsteps, years, spinup_years, nspin))
    log_line(@sprintf("  Julia threads %d   FFT planning threads (JULIA_CPU_THREADS) %d", Threads.nthreads(), Sys.CPU_THREADS))
    log_line(@sprintf("  progress every %.0f model days; output → %s", progress_every_days, file))
    log_line("="^78)
    ndiag = max(1, round(Int, diag_every_days * DAY / dt))
    nspec = max(1, round(Int, spectra_every_days * DAY / dt))
    nsnap = isfinite(snapshot_every_days) ? max(1, round(Int, snapshot_every_days * DAY / dt)) : typemax(Int)
    nvort = isfinite(vorticity_every_days) ? max(1, round(Int, vorticity_every_days * DAY / dt)) : typemax(Int)
    vort1 = Array{Float32, 2}[]
    psi2 = Array{Float32, 2}[]        # lower-layer streamfunction (eddy part), same times
    vorttimes = Float64[]

    ts = Dict{String, Vector{Float64}}()
    push_ts!(k, v) = push!(get!(ts, k, Float64[]), Float64(v))
    spec_sum, k_spec, nspec_acc = nothing, Float64[], 0
    flux_sum = zeros(2)
    snapshots = Array{Float64, 3}[]
    snaptimes = Float64[]
    closure_snapshots = Array{Float64, 2}[]   # eddy-energy field of the closure at the snapshot times
    Fqh = similar(prob.sol)

    t0 = time()
    stopped_early = false
    for n in 1:nsteps
        stepforward!(prob)
        n == nspin && log_line(@sprintf("[%s] spin-up finished at t = %.2f yr; time averaging starts", name, prob.clock.t / YEAR))
        need_fields = has_prognostic_state(cl) || n % ndiag == 0 || n % nvort == 0 || n % nprog == 0 ||
                      (n > nspin && (n % nspec == 0 || n % nsnap == 0))
        need_fields || continue
        f = layerfields(prob.sol, m)
        step_closure!(cl, f, m, dt)

        if n % nprog == 0
            if any(isnan, prob.sol)
                log_line(@sprintf("[%s] NaN detected at step %d (t = %.2f yr), stopping", name, n, prob.clock.t / YEAR))
                stopped_early = true
                break
            end
            log_progress(name, n, nsteps, nspin, prob, f, m, dt, t0)
        end

        if n % ndiag == 0
            if any(isnan, prob.sol)
                log_line(@sprintf("[%s] NaN detected at step %d (t = %.2f yr), stopping", name, n, prob.clock.t / YEAR))
                stopped_early = true
                break
            end
            e = energy_summary(prob, m)
            push_ts!("t", prob.clock.t)
            push_ts!("KE1", e.KE1); push_ts!("KE2", e.KE2); push_ts!("PE", e.PE)
            if !(cl isa NoClosure)
                pv_tendency!(Fqh, cl, f, m)
                push_ts!("closure_energy_rate", closure_energy_rate(Fqh, f, m))
            end
            for (k, v) in pairs(closure_diagnostics(cl, m))
                push_ts!(String(k), v)
            end
            vq = pv_flux(f)
            pf = param_pv_flux(cl, f, m)
            push_ts!("vq_resolved_1", vq[1]); push_ts!("vq_resolved_2", vq[2])
            push_ts!("vq_param_1", pf[1]);    push_ts!("vq_param_2", pf[2])
        end

        if n % nvort == 0
            push!(vort1, Float32.(Array(view(f.ζ, :, :, 1))))
            push!(psi2, Float32.(Array(view(f.ψ, :, :, 2))))
            push!(vorttimes, prob.clock.t)
        end

        if n > nspin && n % nspec == 0
            k_spec, sp = ke_spectra(f, m)
            spec_sum = isnothing(spec_sum) ? sp : spec_sum .+ sp
            flux_sum .+= pv_flux(f)
            nspec_acc += 1
        end

        if n > nspin && n % nsnap == 0
            push!(snapshots, Array(f.q))
            push!(snaptimes, prob.clock.t)
            s = closure_state(cl)
            s === nothing || push!(closure_snapshots, Float64.(s))
        end

    end

    log_line(@sprintf("[%s] %s after %s wall time; writing %s", name,
                      stopped_early ? "stopped early" : "finished", fmt_hms(time() - t0), file))
    jldsave(file;
            params = pp, nx, dt, closure = String(closure),
            closure_settings = closure_settings(cl), timeseries = ts,
            k = k_spec,
            ke_spectrum_mean = isnothing(spec_sum) ? nothing : spec_sum ./ nspec_acc,
            pv_flux_mean = flux_sum ./ max(nspec_acc, 1),
            q_final = Array(layerfields(prob.sol, m).q),
            snapshots, snaptimes, closure_snapshots,
            vorticity_upper = vort1, vorticity_times = vorttimes, psi_lower = psi2,
            metadata = Dict{String, Any}(metadata), spinup_years, aliased_fraction,
            stepper, stopped_early)
    log_line(@sprintf("[%s] saved", name))
    return file
end

# ----------------------------------------------------------------------------
# Progress output. Every line is flushed, so logs of batch jobs (stdout
# redirected to a file, which Julia buffers) stay up to date.
# ----------------------------------------------------------------------------

log_line(s) = (println(stdout, s); flush(stdout); nothing)

function fmt_hms(s)
    s = round(Int, s)
    return @sprintf("%d:%02d:%02d", s ÷ 3600, (s % 3600) ÷ 60, s % 60)
end

"""
One status line: step, progress, model time and phase, time step, advective
CFL number max(|u|+|v|)·dt/Δx (imposed flow included), energies, speed and
estimated remaining wall time.
"""
function log_progress(name, n, nsteps, nspin, prob, f, m, dt, t0)
    umax = maximum(maximum(abs.(view(f.u, :, :, j) .+ m.U[j])) for j in 1:2)
    vmax = maximum(abs.(f.v))
    cfl = (umax + vmax) * dt / m.Δx
    e = energy_summary(prob, m)
    elapsed = time() - t0
    rate = n / max(elapsed, 1e-9)
    eta = (nsteps - n) / rate
    phase = n <= nspin ? "spin-up" : "averaging"
    log_line(@sprintf("[%s] step %d/%d (%5.1f%%) | t = %6.2f yr (%s) | dt = %.0f s | CFL = %.3f | KE1 = %.3e KE2 = %.3e PE = %.3e | %.1f steps/s | elapsed %s | ETA %s",
                      name, n, nsteps, 100 * n / nsteps, prob.clock.t / YEAR, phase, dt, cfl,
                      e.KE1, e.KE2, e.PE, rate, fmt_hms(elapsed), fmt_hms(eta)))
end
