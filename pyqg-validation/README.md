# pyqg validation of the Perezhogin & Glazunov (2023) Reynolds-stress closure

Reproduces Fig. 10 of Perezhogin & Glazunov (2023, JAMES, doi:10.1029/2023MS003771) in pyqg itself,
with the authors' own parameterization code. The diagnostics are identical to the QG testbed
(`QGTestbed/scripts/perezhogin_validation_spectra.jl`), so the two can be overlaid. This separates
differences in the closure from differences in the base model (pyqg vs GeophysicalFlows.jl).

## Contents

| File | Purpose |
|---|---|
| `parameterizations.py` | Perezhogin's modified pyqg parameterization file (`Reynolds_stress`), from the Zenodo record [10.5281/zenodo.8277231](https://zenodo.org/records/8277231), `pyqg.tar.gz`, unchanged |
| `run_pyqg.py` | one run: `--member truth` (256²), `none` (64²) or `reynolds` (64², `Reynolds_stress(FGR=2, Csim=7)`) |
| `analyze_pyqg.py` | spectra, energy transfers, subgrid transfers (SGS Gauss / Sharp), CSV files and the Fig. 10 figure |
| `pyqg_validation.sh` | Slurm job: all 30 runs in parallel, then the analysis |

## Setup (once)

```bash
pip install pyqg matplotlib        # numpy comes with pyqg; pyfftw is optional (faster FFTs)
```

## Run

```bash
sbatch pyqg_validation.sh                  # from this folder
```

Settings are at the top of the script: 10 members, 20 years, statistics over years 5–20, snapshots
every 30 days, C_R = 7, 256² truth and 64² coarse grids. The configuration is pyqg's "eddy" case
(Ross et al. 2023, Table 1): L = 1000 km, β = 1.5e−11 m⁻¹s⁻¹, r_d = 15 km, δ = 0.25, H₁ = 500 m,
U₁ = 0.025 m/s, U₂ = 0, r_ek = 5.787e−7 s⁻¹, Δt = 1 h. pyqg's own numerics are used (AB3, exponential
filter, no dealiasing).

Single run or analysis by hand:

```bash
python run_pyqg.py --member reynolds --seed 1 --years 20 --spinup 5
python analyze_pyqg.py --dir results
python analyze_pyqg.py --dir results --testbed ../QGTestbed/perezhogin-validation/aliased_fraction_0
```

## Output (in `results/`)

- `49N_nx<nx>_<member>_seed<seed>.npz`: PV snapshots after the spin-up and the model parameters
- `spectra_<member>.csv`: k, E, T_KE, T_APE (and P_sgs for `reynolds`)
- `sgs_truth.csv`: subgrid transfer of the coarse-grained truth, Gaussian and sharp operators
- `pyqg_validation.png`: the figure; with `--testbed` the testbed curves are dotted

## Comparing with the testbed

The testbed's validation job (`QGTestbed/scripts/perezhogin_validation.sh`) writes the same CSV files
next to its figure. Pass that folder with `--testbed` (or set `TESTBED=` in the job script).

- **pyqg matches the paper, testbed matches pyqg:** the testbed is validated.
- **pyqg matches the paper, testbed does not:** the difference is in the base numerics of the
  testbed (dealiasing, time stepping, filter), not in the closure; the closure code was checked line by
  line against `Reynolds_stress` and is identical.
- **pyqg does not match the paper:** check the run length, ensemble size and diagnostics first.

Spectra are normalized so that their sum times dk is the domain mean (k is the angular wavenumber).
If the paper used pyqg's own spectral normalization, absolute values may differ from its figure by a
constant factor; shapes and the ratios between curves are comparable.
