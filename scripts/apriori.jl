# Offline (a-priori) evaluation on truth snapshots:
#   julia --project=. scripts/apriori.jl output/49N_nx256_truth_seed1.jld2 64
#
# For each saved truth PV snapshot: coarse-grain to nx_coarse, diagnose the
# exact subgrid PV forcing, and compare with the stateless closures
# (Reynolds backscatter, NEMO default). Prints correlation, relative RMSE and
# the energy transfer of truth vs model.
using QGTestbed
using QGTestbed: physical_parameters, build_problem, make_closure, apriori_skill
using GeophysicalFlows, FourierFlows, JLD2, Statistics, Printf

file, nx_c = ARGS[1], parse(Int, ARGS[2])
d = load(file)
pp, nx_t = d["params"], d["nx"]

prob_hi, m_hi, _ = build_problem(pp, nx_t)
prob_lo, m_lo, _ = build_problem(pp, nx_c)

for name in (:reynolds, :nemo)
    cl = make_closure(name)
    QGTestbed.init_closure!(cl, m_lo)
    res = []
    for q in d["snapshots"]
        MultiLayerQG.set_q!(prob_hi, q)
        push!(res, apriori_skill(cl, prob_hi, m_hi, m_lo))
    end
    @printf("%-10s lat %2.0f  nx %3d  dx/Ld %.2f | corr = [%.2f, %.2f]  rel.RMSE = [%.2f, %.2f]  dE/dt true %+.2e  model %+.2e\n",
            name, pp.lat, nx_c, (pp.L / nx_c) / pp.Ld,
            mean(r.corr[1] for r in res), mean(r.corr[2] for r in res),
            mean(r.rel_rmse[1] for r in res), mean(r.rel_rmse[2] for r in res),
            mean(r.energy_transfer_true for r in res), mean(r.energy_transfer_model for r in res))
end
