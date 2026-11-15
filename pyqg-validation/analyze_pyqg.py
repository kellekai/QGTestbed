"""
Spectra, energy transfers and the Fig. 10 figure of Perezhogin & Glazunov (2023) from the pyqg runs
(run_pyqg.py), with exactly the diagnostics of the QG testbed (perezhogin_validation_spectra.jl),
so both can be overlaid.

    python analyze_pyqg.py --dir results
    python analyze_pyqg.py --dir results --testbed ../QGTestbed/perezhogin-validation/aliased_fraction_0

Diagnostics (time mean over all snapshots and seeds; densities per unit wavenumber in shells of
width 2π/L, normalized so that sum × dk = domain mean; k = angular wavenumber, rad/m):
  E(k)      KE spectrum, thickness weighted, Σ_j (H_j/H) ½|k ψ̂_j|²
  T_KE(k)   resolved KE transfer   Σ_j (H_j/H) Re[ψ̂_j* Ĵ(ψ_j, ζ_j)],  Ĵ(ψ, ζ) = ∇·(u ζ)
  T_APE(k)  resolved APE transfer  (H₁F₁/H) Re[(ψ̂₁ − ψ̂₂)* Ĵ(ψ₁, ψ₂)]
  P(k)      transfer of total energy by the parameterization  −Σ_j (H_j/H) Re[ψ̂_j* Ŝ_j]  (reynolds)
  P_gauss   same for the exact subgrid forcing of the coarse-grained truth, Gaussian exp(−k²(2Δx)²/24)
            + spectral truncation ("SGS Gauss", Ross et al. operator 2)
  P_sharp   same with pyqg's filter exp(−23.6(κ* − 0.65π)⁴) + truncation ("SGS Sharp", operator 1)
No dealiasing mask is applied (pyqg has none).

Writes <dir>/spectra_<member>.csv, <dir>/sgs_truth.csv and <dir>/pyqg_validation.png. With
--testbed DIR, the testbed's CSV files (written by plot_perezhogin_validation.jl) are overlaid.
"""
import argparse
import csv
import glob
import os
import re
import numpy as np

NAME = re.compile(r"^(\d+)N_nx(\d+)_(truth|none|reynolds)_seed(\d+)\.npz$")


# --- spectral helpers (numpy rfft2 over the last two axes (y, x)) ---------------------------
class Grid:
    def __init__(self, n, L):
        self.n, self.L = n, L
        self.k = 2 * np.pi / L * np.fft.rfftfreq(n, 1.0 / n)[None, :]     # x, last axis
        self.l = 2 * np.pi / L * np.fft.fftfreq(n, 1.0 / n)[:, None]      # y
        self.K2 = self.k**2 + self.l**2
        self.w = np.full(self.K2.shape, 2.0)
        self.w[:, 0] = 1.0
        if n % 2 == 0:
            self.w[:, -1] = 1.0
        self.dk = 2 * np.pi / L
        self.dx = L / n

    fft = staticmethod(lambda a: np.fft.rfft2(a))

    def ifft(self, ah):
        return np.fft.irfft2(ah, s=(self.n, self.n))

    def cross(self, ah, bh):
        """Per-mode density of the domain mean of a·b."""
        return self.w * np.real(np.conj(ah) * bh) / self.n**4

    def shells(self, D):
        K = np.sqrt(self.K2)
        nb = self.n // 2
        b = np.rint(K / self.dk).astype(int)
        ok = b <= nb
        s = np.bincount(b[ok].ravel(), weights=D[ok].ravel(), minlength=nb + 1)[: nb + 1]
        return s / self.dk

    def wavenumbers(self):
        return np.arange(self.n // 2 + 1) * self.dk


def invert(qh, g, F1, F2):
    """ψ̂ from q̂ for q₁ = ∇²ψ₁ + F₁(ψ₂ − ψ₁), q₂ = ∇²ψ₂ + F₂(ψ₁ − ψ₂)."""
    a11, a12, a21, a22 = -g.K2 - F1, F1, F2, -g.K2 - F2
    det = a11 * a22 - a12 * a21
    det = np.where(det == 0, np.inf, det)
    p1 = (a22 * qh[0] - a12 * qh[1]) / det
    p2 = (-a21 * qh[0] + a11 * qh[1]) / det
    return np.array([p1, p2])


def fields(qh, g, F1, F2):
    ph = invert(qh, g, F1, F2)
    zh = -g.K2 * ph
    u = g.ifft(-1j * g.l * ph)
    v = g.ifft(1j * g.k * ph)
    return dict(qh=qh, ph=ph, zh=zh, q=g.ifft(qh), u=u, v=v, z=g.ifft(zh), p=g.ifft(ph))


def div_flux(fx, fy, g):
    """Spectral ∇·(fx, fy) of physical fields."""
    return 1j * g.k * g.fft(fx) + 1j * g.l * g.fft(fy)


def resolved_spectra(f, g, H, F1):
    Hs = sum(H)
    E = np.zeros(g.K2.shape); TKE = np.zeros(g.K2.shape)
    for j in range(2):
        J = div_flux(f["u"][j] * f["z"][j], f["v"][j] * f["z"][j], g)
        E += H[j] / Hs * (-0.5) * g.cross(f["ph"][j], f["zh"][j])
        TKE += H[j] / Hs * g.cross(f["ph"][j], J)
    J12 = div_flux(f["u"][0] * f["p"][1], f["v"][0] * f["p"][1], g)
    TAPE = H[0] * F1 / Hs * g.cross(f["ph"][0] - f["ph"][1], J12)
    return g.shells(E), g.shells(TKE), g.shells(TAPE)


def forcing_transfer(f, Sh, g, H):
    Hs = sum(H)
    P = np.zeros(g.K2.shape)
    for j in range(2):
        P -= H[j] / Hs * g.cross(f["ph"][j], Sh[j])
    return g.shells(P)


def truncate(Ah, ghi, glo):
    """Spectral truncation of a truth-grid rfft2 field (..., n_hi, n_hi/2+1) to the coarse grid."""
    nlo, nhi = glo.n, ghi.n
    rows = np.r_[0:nlo // 2, nhi - nlo // 2:nhi]
    B = Ah[..., rows, : nlo // 2 + 1] * (nlo / nhi) ** 2
    B[..., :, -1] = 0          # coarse Nyquist column
    B[..., nlo // 2, :] = 0    # coarse Nyquist row
    return B


def sgs_forcing(qh_hi, ghi, glo, G, F1, F2):
    """Coarse-grained fields and the exact subgrid PV forcing S = J(ψ̄, q̄) − \\overline{J(ψ, q)}."""
    fh = fields(qh_hi, ghi, F1, F2)
    Jh = np.array([div_flux(fh["u"][j] * fh["q"][j], fh["v"][j] * fh["q"][j], ghi) for j in range(2)])
    Jbar = truncate(G * Jh, ghi, glo)
    qbarh = truncate(G * qh_hi, ghi, glo)
    fl = fields(qbarh, glo, F1, F2)
    Jres = np.array([div_flux(fl["u"][j] * fl["q"][j], fl["v"][j] * fl["q"][j], glo) for j in range(2)])
    return fl, Jres - Jbar


def sharp_transfer(ghi, dx_lo):
    kap = np.sqrt(ghi.K2) * dx_lo
    kc = 0.65 * np.pi
    return np.where(kap < kc, 1.0, np.exp(-23.6 * (kap - kc) ** 4))


# --- analysis -------------------------------------------------------------------------------
def analyze(folder, nx_coarse):
    groups = {}
    for f in sorted(glob.glob(os.path.join(folder, "*.npz"))):
        m = NAME.match(os.path.basename(f))
        if m:
            groups.setdefault(m.group(3), []).append(f)
    if "truth" not in groups:
        raise SystemExit(f"no truth runs in {folder}")
    res, meta = {}, {}
    for member in ("truth", "none", "reynolds"):
        if member not in groups:
            print(f"no '{member}' runs, skipped")
            continue
        acc, nsnap = {}, 0
        add = lambda k, v: acc.__setitem__(k, acc[k] + v if k in acc else v.copy())
        for file in groups[member]:
            d = np.load(file)
            n, L = int(d["nx"]), float(d["L"])
            H, F1, F2 = (float(d["H1"]), float(d["H2"])), float(d["F1"]), float(d["F2"])
            meta.update(L=L, rd=float(d["rd"]), nx_coarse=nx_coarse)
            g = Grid(n, L)
            if member == "reynolds":
                import sys
                sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
                from parameterizations import Reynolds_stress
                param = Reynolds_stress(FGR=float(d["fgr"]), Csim=float(d["csim"]))
                meta["csim"] = float(d["csim"])
            if member == "truth":
                glo = Grid(nx_coarse, L)
                G_gauss = np.exp(-g.K2 * (2 * glo.dx) ** 2 / 24)
                G_sharp = sharp_transfer(g, glo.dx)
            for q in d["q"]:
                qh = g.fft(q.astype(np.float64))
                f = fields(qh, g, F1, F2)
                E, TKE, TAPE = resolved_spectra(f, g, H, F1)
                add("E", E); add("T_KE", TKE); add("T_APE", TAPE)
                if member == "reynolds":
                    S = param(PyqgLike(f, g, H))
                    add("P_sgs", forcing_transfer(f, g.fft(S), g, H))
                elif member == "truth":
                    fl, S = sgs_forcing(qh, g, glo, G_gauss, F1, F2)
                    add("P_gauss", forcing_transfer(fl, S, glo, H))
                    fl, S = sgs_forcing(qh, g, glo, G_sharp, F1, F2)
                    add("P_sharp", forcing_transfer(fl, S, glo, H))
                nsnap += 1
            print(f"{os.path.basename(file):36s} {len(d['q']):4d} snapshots")
        out = {k: v / nsnap for k, v in acc.items()}
        out.update(k=g.wavenumbers(), nx=n, nsnap=nsnap, nseeds=len(groups[member]))
        if member == "truth":
            out["k_coarse"] = glo.wavenumbers()
        res[member] = out
        print(f"{member:8s} {n}²  KE = {out['E'].sum() * g.dk:.4e} m²/s²   ΣT_KE·dk = {out['T_KE'].sum() * g.dk:+.2e}"
              f"   ΣT_APE·dk = {out['T_APE'].sum() * g.dk:+.2e} m²/s³")
    return res, meta


class PyqgLike:
    """The attributes of a pyqg model that Reynolds_stress uses, built from a snapshot, so the
    parameterization is evaluated with Perezhogin's own code. pyqg arrays are (nz, ny, nx),
    spectral (nz, ny, nx/2+1), with k along the last axis, as here."""
    def __init__(self, f, g, H):
        self.q, self.u, self.v = f["q"], f["u"], f["v"]
        self.k = np.broadcast_to(g.k, g.K2.shape)
        self.l = np.broadcast_to(g.l, g.K2.shape)
        self.dx = g.dx
        self.Hi = np.array(H)
        self._g = g

    def fft(self, a):
        return self._g.fft(a)

    def ifft(self, ah):
        return self._g.ifft(ah)


def write_csv(path, cols):
    names = list(cols)
    n = len(cols[names[0]])
    with open(path, "w", newline="") as io:
        w = csv.writer(io)
        w.writerow(names)
        for i in range(n):
            w.writerow([f"{cols[c][i]:.8e}" for c in names])


def read_csv(path):
    with open(path) as io:
        r = csv.DictReader(io)
        rows = list(r)
    return {c: np.array([float(x[c]) for x in rows]) for c in rows[0]} if rows else {}


def plot(res, meta, out_png, testbed=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    nxc = meta["nx_coarse"]
    kmax = np.pi * nxc / meta["L"]
    xl = (1e-5, 1.1 * kmax)
    kLd = 1 / meta["rd"]
    nt = res["truth"]["nx"]
    style = {"none": dict(color="0.5", lw=2, ls="-", label=f"{nxc}²"),
             "reynolds": dict(color="purple", lw=2.5, ls="--", label=f"{nxc}², Reynolds"),
             "truth": dict(color="k", lw=1, ls="-", label=f"{nt}²")}
    order = [m for m in ("none", "reynolds", "truth") if m in res]
    fig, axs = plt.subplots(2, 2, figsize=(11, 8.5))
    (a, b), (c, d) = axs

    def sel(k):
        return (k > 0) & (k <= xl[1])

    for ax, key, title in ((a, "E", "(a) KE spectrum"), (b, "T_KE", "(b) Resolved KE transfer"),
                           (c, "T_APE", "(c) Resolved APE transfer")):
        scale = 1.0 if key == "E" else 1e6
        for mbr in order:
            r = res[mbr]; i = sel(r["k"])
            ax.plot(r["k"][i], r[key][i] * scale, **style[mbr])
            if testbed and mbr in testbed and key in testbed[mbr]:
                t = testbed[mbr]; j = sel(t["k"])
                ax.plot(t["k"][j], t[key][j] * scale, color=style[mbr]["color"], lw=1.2, ls=":",
                        label=style[mbr]["label"] + " (testbed)")
        ax.set_title(title)
        ax.axvline(kLd, color="steelblue", alpha=0.6)
        if key == "E":
            ax.set_yscale("log"); ax.set_ylabel("E(k), m³/s²")
            top = max(np.nanmax(res[m]["E"][sel(res[m]["k"])]) for m in order)
            ax.set_ylim(top / 150, top * 2)
            k0, k1 = 8e-5, 1.5e-4
            ax.plot([k0, k1], [0.8 * top, 0.8 * top * (k1 / k0) ** -3], "k-", lw=1)
            ax.text(1.25e-4, 0.75 * top, "k⁻³")
            ax.legend(loc="lower left", fontsize=8)
        else:
            ax.axhline(0, color="lightgray", ls=":")
            ax.set_ylabel("×10⁻⁶ m³/s³")
    d.set_title("(d) Subgrid KE+APE transfer")
    d.axhline(0, color="lightgray", ls=":")
    d.axvline(kLd, color="steelblue", alpha=0.6)
    if "reynolds" in res:
        r = res["reynolds"]; i = sel(r["k"])
        d.plot(r["k"][i], r["P_sgs"][i] * 1e6, **dict(style["reynolds"], label=f"{nxc}², Reynolds (C_R = {meta.get('csim', 7):g})"))
        if testbed and "reynolds" in testbed and "P_sgs" in testbed["reynolds"]:
            t = testbed["reynolds"]; j = sel(t["k"])
            d.plot(t["k"][j], t["P_sgs"][j] * 1e6, color="purple", lw=1.2, ls=":", label="Reynolds (testbed)")
    t = res["truth"]; kc = t["k_coarse"]; i = kc > 0
    d.plot(kc[i], t["P_gauss"][i] * 1e6, "k-.", lw=1.3, label=f"SGS Gauss ({nt}² → {nxc}²)")
    d.plot(kc[i], t["P_sharp"][i] * 1e6, "k--", lw=1.3, label=f"SGS Sharp ({nt}² → {nxc}²)")
    d.set_ylabel("×10⁻⁶ m³/s³")
    d.legend(loc="lower left", fontsize=8)
    for ax in axs.flat:
        ax.set_xscale("log"); ax.set_xlim(*xl); ax.set_xlabel("wavenumber, m⁻¹")
    if testbed:
        b.legend(loc="upper right", fontsize=7)
    fig.suptitle(f"pyqg with Perezhogin's Reynolds_stress: eddy configuration, rd = {meta['rd'] / 1e3:.0f} km"
                 + ("  (dotted: QG testbed)" if testbed else ""))
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    print(f"saved {out_png}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", default="results")
    p.add_argument("--nx-coarse", type=int, default=64)
    p.add_argument("--testbed", help="folder with the testbed's spectra_<member>.csv and sgs_truth.csv")
    a = p.parse_args()

    res, meta = analyze(a.dir, a.nx_coarse)
    for mbr, r in res.items():
        cols = {k: r[k] for k in ("k", "E", "T_KE", "T_APE", "P_sgs") if k in r}
        write_csv(os.path.join(a.dir, f"spectra_{mbr}.csv"), cols)
    t = res["truth"]
    write_csv(os.path.join(a.dir, "sgs_truth.csv"), {"k_coarse": t["k_coarse"], "P_gauss": t["P_gauss"], "P_sharp": t["P_sharp"]})

    testbed = None
    if a.testbed:
        testbed = {}
        for mbr in ("truth", "none", "reynolds"):
            f = os.path.join(a.testbed, f"spectra_{mbr}.csv")
            if os.path.isfile(f):
                testbed[mbr] = read_csv(f)
        if not testbed:
            print(f"no testbed CSV files in {a.testbed}; run plot_perezhogin_validation.jl first")
            testbed = None
    plot(res, meta, os.path.join(a.dir, "pyqg_validation.png"), testbed)


if __name__ == "__main__":
    main()
