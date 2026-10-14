# Quick consistency check, run this first:   julia --project=. scripts/smoke_test.jl
#
# 1. Builds every closure at 49°N on a 64² grid and runs 1 model year.
# 2. Checks: no NaNs, energy grows from the random initial state, and the sign
#    of each closure's energy transfer (NEMO, GEOMETRIC: < 0; Reynolds
#    backscatter: > 0 most of the time; MEKE: mixed, GM/biharmonic vs backscatter).
# 3. Checks that the GM tendency removes PE at exactly the rate the eddy-energy
#    budgets receive (energy consistency of the coupling).
using QGTestbed
using QGTestbed: physical_parameters, CaseConfig, build_problem, make_closure,
                 set_random_initial_condition!, layerfields, pv_tendency!,
                 closure_energy_rate, interface_gradient, gm_pv_tendency!, gm_energy_release,
                 energy_summary
using GeophysicalFlows, FourierFlows, Statistics, Printf

latitude_table()
println()

pp = physical_parameters(CaseConfig(lat = 49.0))

# --- 3. GM energy consistency on a spun-up-ish random state -------------------
prob, m, _ = build_problem(pp, 64)
set_random_initial_condition!(prob, pp; amplitude = 1e-6)
for _ in 1:200; stepforward!(prob); end
f = layerfields(prob.sol, m)
# With constant κ the resolved eddy PE loss must equal the eddy part of the
# release term exactly (to round-off):
κ0 = fill!(similar(m.real2), 500.0)
Fqh0 = zero(prob.sol)
gxa, gya = interface_gradient(f.ψh[:, :, 1] .- f.ψh[:, :, 2], m)
gm_pv_tendency!(Fqh0, κ0, gxa, gya, m)
gx_eddy = gxa; gy_eddy = gya .+ m.ΔU
@printf("GM energy check (constant κ): resolved dE/dt = %.4e, -release/H = %.4e\n",
        closure_energy_rate(Fqh0, f, m), -mean(gm_energy_release(κ0, gx_eddy, gy_eddy, m)) / m.Htot)

# --- 1–2. one year per closure --------------------------------------------------
for name in (:none, :nemo, :geometric, :geometric_nosplit, :meke, :reynolds)
    local prob, m, f
    cl = make_closure(name)
    prob, m, cl = build_problem(pp, 64; closure = cl)
    set_random_initial_condition!(prob, pp)
    E0 = energy_summary(prob, m).Etot
    nsteps = round(Int, 365 * 86400 / prob.clock.dt)
    rates = Float64[]
    ok = true
    for n in 1:nsteps
        stepforward!(prob)
        f = layerfields(prob.sol, m)
        QGTestbed.step_closure!(cl, f, m, prob.clock.dt)
        if n % 24 == 0
            if any(isnan, prob.sol)
                ok = false; break
            end
            if !(cl isa NoClosure)
                Fh = similar(prob.sol)
                pv_tendency!(Fh, cl, f, m)
                push!(rates, closure_energy_rate(Fh, f, m))
            end
        end
    end
    E1 = energy_summary(prob, m).Etot
    @printf("%-18s  NaN-free: %-5s  E(0) = %.2e  E(1 yr) = %.2e  mean closure dE/dt = %+.2e  diagnostics = %s\n",
            name, ok, E0, E1, isempty(rates) ? 0.0 : mean(rates),
            string(QGTestbed.closure_diagnostics(cl, m)))
end
