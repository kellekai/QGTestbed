"""
One pyqg run of the Perezhogin & Glazunov (2023) Fig. 10 case: the pyqg "eddy" configuration
(Ross et al. 2023, Table 1), with their own Reynolds-stress parameterization
(parameterizations.py from the Zenodo record 10.5281/zenodo.8277231, same folder).

    python run_pyqg.py --member truth    --seed 1      # 256², no parameterization
    python run_pyqg.py --member none     --seed 1      # 64², no parameterization
    python run_pyqg.py --member reynolds --seed 1      # 64², Reynolds_stress(FGR=2, Csim=7)

Options: --years 20 --spinup 5 --snapshot-days 30 --nx-truth 256 --nx-coarse 64 --csim 7
         --fgr 2 --outdir results --threads 1

Writes <outdir>/49N_nx<nx>_<member>_seed<seed>.npz with the PV snapshots q (float32,
(nsnap, 2, ny, nx)) after the spin-up and the model parameters needed by analyze_pyqg.py.
Uses pyqg's own numerics (AB3, exponential filter, no dealiasing), dt = 1 h as in Ross et al.
"""
import argparse
import os
import time
import numpy as np
import pyqg

YEAR = 365 * 86400.0
DAY = 86400.0

# pyqg "eddy" configuration (Ross et al. 2023, Table 1; these are also pyqg's QGModel defaults)
EDDY = dict(L=1e6, beta=1.5e-11, rek=5.787e-7, rd=15000.0, delta=0.25, H1=500.0, U1=0.025, U2=0.0)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--member", required=True, choices=["truth", "none", "reynolds"])
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--years", type=float, default=20.0)
    p.add_argument("--spinup", type=float, default=5.0)
    p.add_argument("--snapshot-days", type=float, default=30.0)
    p.add_argument("--nx-truth", type=int, default=256)
    p.add_argument("--nx-coarse", type=int, default=64)
    p.add_argument("--csim", type=float, default=7.0, help="C_R of the Reynolds parameterization")
    p.add_argument("--fgr", type=float, default=2.0, help="filter-to-grid ratio of the Reynolds parameterization")
    p.add_argument("--dt", type=float, default=3600.0)
    p.add_argument("--threads", type=int, default=1, help="FFT threads (pyfftw), pyqg's ntd")
    p.add_argument("--outdir", default="results")
    a = p.parse_args()

    nx = a.nx_truth if a.member == "truth" else a.nx_coarse
    param = None
    if a.member == "reynolds":
        from parameterizations import Reynolds_stress      # Perezhogin's file, same folder
        param = Reynolds_stress(FGR=a.fgr, Csim=a.csim)

    m = pyqg.QGModel(nx=nx, dt=a.dt, tmax=a.years * YEAR, twrite=int(10 * YEAR / a.dt),
                     tavestart=a.years * YEAR, ntd=a.threads, q_parameterization=param, **EDDY)

    # reproducible small random initial PV, as in the testbed
    rng = np.random.default_rng(a.seed)
    q0 = 1e-7 * rng.standard_normal((2, nx, nx))
    q0 -= q0.mean(axis=(1, 2), keepdims=True)
    m.set_q(q0)

    os.makedirs(a.outdir, exist_ok=True)
    name = f"49N_nx{nx}_{a.member}_seed{a.seed}"
    print(f"[{name}] pyqg {pyqg.__version__}, nx = {nx}, dt = {a.dt:.0f} s, {a.years:g} yr, "
          f"snapshots every {a.snapshot_days:g} d after {a.spinup:g} yr, parameterization: {param}", flush=True)

    snaps, times = [], []
    t0 = time.time()
    tsnapint = a.snapshot_days * DAY
    for _ in m.run_with_snapshots(tsnapstart=a.spinup * YEAR, tsnapint=tsnapint):
        snaps.append(m.q.astype(np.float32).copy())
        times.append(m.t)
        if len(snaps) % 12 == 1:
            ke = 0.5 * np.mean(m.u**2 + m.v**2)
            el = time.time() - t0
            print(f"[{name}] t = {m.t / YEAR:6.2f} yr  snapshots {len(snaps):4d}  KE = {ke:.3e} m²/s²  "
                  f"elapsed {el / 60:6.1f} min", flush=True)

    file = os.path.join(a.outdir, name + ".npz")
    np.savez_compressed(file, q=np.array(snaps), t=np.array(times), nx=nx, L=m.L,
                        H1=m.Hi[0], H2=m.Hi[1], F1=m.F1, F2=m.F2, beta=m.beta, rd=m.rd, delta=m.delta,
                        U1=EDDY["U1"], U2=EDDY["U2"], rek=m.rek, dt=a.dt, member=a.member, seed=a.seed,
                        csim=a.csim if param is not None else np.nan, fgr=a.fgr, years=a.years,
                        spinup=a.spinup, pyqg_version=pyqg.__version__)
    print(f"[{name}] saved {file} ({len(snaps)} snapshots, {(time.time() - t0) / 60:.1f} min)", flush=True)


if __name__ == "__main__":
    main()
