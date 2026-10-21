# QGTestbed — comparing mesoscale eddy parametrizations across latitudes

This repository compares four mesoscale eddy parametrizations, **as published and
without retuning**, in one idealized model. The question is how each behaves as
the ratio of grid spacing to deformation radius (Δx/Ld) changes with latitude.

| key | parametrization | scale aware |
|---|---|---|
| `nemo` | NEMO 4 default: bilaplacian viscosity + constant GM | no |
| `reynolds` | Perezhogin & Glazunov (2023): Reynolds-stress backscatter (Germano decomposition) | no (fixed C_R) |
| `geometric` | Mak et al. (2023): GEOMETRIC with the "splitting" filter | yes |
| `meke` | Jansen et al. (2019): MEKE with GM + backscatter, resolution function | yes |
| `none` | reference: no closure | – |
| `geometric_nosplit` | reference: original GEOMETRIC without splitting | no |

Every closure is judged against an eddy-resolving **truth** run at the same latitude.

---

## 1. The model

- 2-layer quasi-geostrophic model, doubly periodic β-plane, 1000 km × 1000 km
  (GeophysicalFlows.jl, `MultiLayerQG`).
- **Forcing:** a fixed eastward flow in the upper layer (U₁ = ΔU, U₂ = 0). This keeps
  the interface tilted (thermal wind) and maintains a background PV gradient that changes
  sign between the layers, which drives baroclinic instability. There is no wind and no
  buoyancy forcing. Energy path: fixed mean APE → eddy APE → eddy KE → bottom drag
  (large scales); enstrophy → grid scale (spectral filter or closure).
- **Latitude sweep:** one run per latitude (20, 30, 40, 50, 60, 70°N) with that latitude's
  f₀ and β. Stratification is fixed, so Ld ∝ 1/f₀. ΔU is scaled to keep the same
  supercriticality ξ = ΔU F₂/β = 1.48 everywhere. The 49°N equivalent is the pyqg "eddy"
  configuration (Ld = 15 km, β = 1.5×10⁻¹¹ m⁻¹s⁻¹), on which Perezhogin & Glazunov tuned C_R = 7.

| lat | Ld [km] | ΔU [m/s] | truth grid | Δx/Ld at 32² | 64² | 128² |
|----:|--------:|---------:|-----------:|------:|-----:|-----:|
| 20 | 33.1 | 0.174 | 256² | 0.94 | 0.47 | 0.24 |
| 30 | 22.6 | 0.075 | 256² | 1.38 | 0.69 | 0.35 |
| 40 | 17.6 | 0.040 | 256² | 1.77 | 0.89 | 0.44 |
| 50 | 14.8 | 0.024 | 512² | 2.11 | 1.06 | 0.53 |
| 60 | 13.1 | 0.014 | 512² | 2.39 | 1.20 | 0.60 |
| 70 | 12.0 | 0.008 | 512² | 2.59 | 1.30 | 0.65 |

Layer depths 500 m / 2000 m, linear bottom drag 5.787×10⁻⁷ s⁻¹.

---

## 2. Installation (once per machine)

Requirements: Julia ≥ 1.10 (a recent 1.11.x is recommended) and GeophysicalFlows ≥ 0.16.

```bash
cd QGTestbed
julia setup.jl              # model environment (this folder)
julia setup.jl plotting     # optional: separate plotting environment in ./plotting (Plots.jl, GR)
```

- The model environment contains no plotting packages, so compute nodes never compile them.
- On a cluster, run `setup.jl` on a login node **with the same Julia module the batch jobs use**,
  so packages are compiled once and jobs don't compile in parallel.
- The package is shipped **without `Project.toml`**: `setup.jl` creates it and adds the dependencies.
  Unpacking an updated version over an existing folder therefore never touches the installed
  environment. **After every update, run `julia setup.jl` again** on the login node: it recompiles
  QGTestbed once, before any batch job does. Don't run a partial `Pkg.add` instead (see Troubleshooting).

## 3. Check the installation (≈ 1–2 min)

```bash
julia --project=. scripts/smoke_test.jl
```

Expected output:
- the latitude table above;
- `GM energy check`: two **identical** numbers;
- for every closure: `NaN-free: true`. The mean closure dE/dt should be negative for `nemo`
  and both GEOMETRIC variants and positive (but tiny) for `reynolds`.

After one model year the flow is still in its growth phase, so energies are small. The test
checks stability and signs, not equilibrium.

Optional check of the output code:
```bash
julia --project=. -e 'using QGTestbed; run_simulation(lat = 50, nx = 32, closure = :meke, years = 0.3, spinup_years = 0.1, outdir = "test")'
```

---

## 4. Run the truth simulations

One eddy-resolving run per latitude (no closure, ≥ 4 grid points per Ld), 20 years, of which
the first 5 are spin-up. They save PV snapshots every 30 days after spin-up (for
coarse-graining) and upper-layer vorticity every 30 days from t = 0.

**Slurm (one node, 112 cores):** edit the account and QoS lines, then
```bash
sbatch scripts/slurm/truth.sh            # seed 1;  SEED=2 sbatch scripts/slurm/truth.sh for more
```
To get a first comparison before the full coarse sweep, `scripts/slurm/truth_and_50N.sh` runs the
same truth runs plus the six 64² closure runs at 50°N (seed 1) in one job. Full listings of all
batch scripts are in the Appendix.

**Workstation:**
```bash
for lat in 50 60 70; do
    JULIA_CPU_THREADS=16 nohup julia -t 16 --project=. scripts/run_member.jl $lat 0 truth 20 5 1 > truth_${lat}N.log 2>&1 &
done
for lat in 20 30 40; do
    JULIA_CPU_THREADS=8 nohup julia -t 8 --project=. scripts/run_member.jl $lat 0 truth 20 5 1 > truth_${lat}N.log 2>&1 &
done
```

Arguments of `run_member.jl`: `lat nx closure years spinup_years [seed] [gpu]`; `nx = 0` selects
the truth resolution automatically.

**Threads.** Only the FFTs are multithreaded, and only the 512² runs profit noticeably.
FourierFlows plans FFTs for `Sys.CPU_THREADS` threads, so always set `JULIA_CPU_THREADS`
equal to `-t`. Benchmark once: compare the wall time per 100 model days printed in the logs
for `-t 1`, `4` and `16`.

## 5. Run the coarse simulations (all parametrizations)

324 runs: 6 latitudes × 3 grids (32², 64², 128²) × 6 closures × 3 seeds, single-threaded.

```bash
julia --project=. scripts/sweep.jl | grep -v truth | sed 's/^/JULIA_CPU_THREADS=1 /' > coarse_jobs.txt
sbatch scripts/slurm/coarse.sh                  # GNU parallel, 100 runs at a time
```

If the 50°N 64² seed-1 runs were already done with `truth_and_50N.sh`, remove them from the list
before submitting, otherwise they run twice:
```bash
grep -v -E 'run_member.jl 50 64 [a-z_]+ 20 5 1( |$)' coarse_jobs.txt > tmp && mv tmp coarse_jobs.txt
```

Each run writes its file only when it finishes. If the job reaches the time limit, submit it
again: `--resume-failed` skips runs that completed successfully and reruns failed ones. Per-run logs are in `logs/`.

To run a single case: `julia --project=. scripts/run_member.jl 50 64 reynolds 20 5 1`.

---

## 6. Output files

`output/<lat>N_nx<n>_<closure>_seed<s>.jld2`, e.g. `50N_nx64_none_seed1.jld2`. Load with
`using JLD2; d = load(file)`.

| key | content |
|---|---|
| `params`, `nx`, `dt`, `closure`, `closure_settings` | setup of the run |
| `timeseries` | daily: `t`, `KE1`, `KE2`, `PE`, `vq_resolved_1/2`, `vq_param_1/2`; closure runs also `closure_energy_rate` and closure diagnostics (`kappa_gm`, `E_param`, `MEKE`, `C_R`, …) |
| `k`, `ke_spectrum_mean` | KE spectrum per layer, time mean after spin-up |
| `pv_flux_mean` | resolved meridional PV flux per layer, time mean after spin-up |
| `vorticity_upper`, `vorticity_times` | upper-layer relative vorticity every 30 days from t = 0 (Float32) |
| `psi_lower` | lower-layer streamfunction at the same times (Float32; runs started after it was added) |
| `snapshots`, `snaptimes` | full PV fields every 30 days after spin-up (truth runs only) |
| `q_final` | final PV field |

Sizes: a few MB for coarse runs, about 500 MB for 512² truth runs (vorticity and lower-layer streamfunction snapshots).

---

## 7. Analysis

All plotting uses the separate environment: `julia --project=plotting …` (run from this folder).

1. **Equilibration and spectra.** Energy time series and mean upper-layer KE spectra, one or many runs overlaid:
   ```bash
   julia --project=plotting scripts/plot_summary.jl output/50N_nx64_none_seed1.jld2
   julia --project=plotting scripts/plot_summary.jl output/50N_nx64_*_seed1.jld2 output/50N_nx512_truth_seed1.jld2 -o output/summary_50N_nx64.png
   ```
   Energies must level off well before year 5, otherwise lengthen the spin-up. The spectra are
   computed from the vorticity snapshots and normalized so that ∫E(k)dk is the domain-mean KE,
   so truth and coarse runs are directly comparable (the stored `ke_spectrum_mean` of runs made
   before this fix used unnormalized FFT amplitudes and should not be compared across grids).
2. **Vorticity panels** (eddy development, truth vs every closure):
   ```bash
   julia --project=plotting scripts/plot_vorticity.jl 50 64              # times chosen automatically
   julia --project=plotting scripts/plot_vorticity.jl 50 64 0.5 1 2 5    # times in years
   julia --project=plotting scripts/plot_vorticity.jl 50 64 -d output/tuned -t output -o figures   # folders
   ```
   Writes `output/vorticity_50N_nx64.png` (`-d` run folder, `-t` truth folder, `-o` figure folder). Columns: truth (coarse-grained to the same grid)
   and each closure; rows: times; one colour scale per row from the truth.
2b. **Mean flow and eddies** (as Jansen et al. 2019, Fig. 3):
   ```bash
   julia --project=plotting scripts/plot_circulation.jl 40 64          # snapshot = last saved time
   julia --project=plotting scripts/plot_circulation.jl 40 64 12.5     # snapshot at 12.5 years
   ```
   Writes `output/circulation_40N_nx64.png`. Contours: sea surface height η = f₀ψ₁/g (the QG
   upper-layer streamfunction, imposed flow included; black > 0, white < 0). Shading: magnitude of the
   vertically integrated transport |H₁u₁ + H₂u₂| on a log scale. Row (a) is the time mean after the
   spin-up, row (b) a snapshot of the deviations from it; columns are the truth (coarse-grained) and each
   closure. In this doubly periodic testbed the mean is mostly the imposed flow, so row (b) is the
   informative one. The transport needs the lower-layer streamfunction (`psi_lower`), saved by runs
   started after this change; for older runs the shading falls back to the upper-layer speed |u₁|.
3. **KE spectra:** compare each closure with the truth spectrum at the same latitude, left of the
   dotted cutoff line in the `plot_summary.jl` figure (the scales the coarse grid resolves).
4. **PV fluxes.** For GM-type closures (`nemo`, `geometric*`, `meke`) compare
   `vq_resolved + vq_param` with the truth's `vq_resolved`. Because the imposed shear is fixed,
   GM never changes the mean state here; its flux is only diagnosed.
5. **Offline (a priori) skill** of the stateless closures on truth snapshots:
   ```bash
   julia --project=. scripts/apriori.jl output/50N_nx512_truth_seed1.jld2 64
   ```
   Prints correlation and relative RMSE between modelled and true subgrid PV forcing, and both
   energy transfers.
6. **Across latitudes (scale awareness).** Reduce every coarse run to three numbers relative to the
   truth at the same latitude and plot them against the resolution parameter Ld/Δx (grid points per
   deformation radius, Hallberg 2013), shaded as unresolved (< 1), gray zone (1–2) and resolved (> 2):
   ```bash
   julia --project=plotting scripts/scale_awareness.jl              # all runs in output/
   julia --project=plotting scripts/scale_awareness.jl -s 5 -d output
   julia --project=plotting scripts/scale_awareness.jl -d output/tuned -t output   # tuned sweep (7c)
   julia --project=plotting scripts/scale_awareness.jl -g 64                       # only the 64² grid
   ```
   - spectral error: RMS of log₁₀(E_model/E_truth) for 0 < k ≤ k_c/2 (k_c = dealiasing cutoff);
     0 is perfect, 1 means a factor 10 off on average, capped at 6 (eddies gone);
   - resolved KE: model KE / truth KE at k ≤ k_c, minus 1;
   - PV flux: (resolved + parametrized upper-layer ⟨vq⟩) / truth ⟨vq⟩, minus 1.

   Writes `output/scale_awareness.csv` (one row per latitude, grid and closure, mean and range over
   seeds) and one figure per metric, each with one row per grid (32², 64², 128²), an Ld/Δx axis fitted to
   each grid and the latitudes on a second axis on top (colours = closures): `output/scale_awareness_spectral_error.png`,
   `output/scale_awareness_resolved_KE.png`, `output/scale_awareness_PV_flux.png`. A scale-aware
   closure stays near 0 across the whole Ld/Δx range; a non-scale-aware one drifts with Ld/Δx.
7. **MEKE energy budget** (Jansen 2019 closure): what sets the subgrid eddy energy e, and with it
   K_GM and the backscatter. Needs `meke` runs made with the current code.
   ```bash
   julia --project=plotting scripts/meke_budget.jl output/40N_nx64_meke_seed1.jld2 output/60N_nx64_meke_seed1.jld2
   julia --project=plotting scripts/meke_budget.jl output/*N_nx64_meke_seed1.jld2 -t output -o output/meke_budget_nx64.png
   ```
   Prints per run the time means after the spin-up of e, √(2e), L_mix, R and K_GM (with the truth's
   κ_eff = −⟨v₁q₁⟩/(F₁ΔU) from the `-t` folder for comparison), the four budget terms of ∂e/∂t
   (GM PE release and biharmonic dissipation as sources, backscatter and drag as sinks) with their
   share of the source, and the drag damping time, including the part set by U_bg alone. The figure
   shows the terms against time. The run files store these terms as `meke_src_GM`,
   `meke_src_biharm`, `meke_sink_bs`, `meke_sink_drag`, `meke_net`, `meke_drag_rate`,
   `meke_drag_rate_Ubg` and `meke_Lmix` in `timeseries`.

8. **Validation against Perezhogin & Glazunov (2023), Fig. 8.** The 49°N member is the pyqg
   "eddy" configuration (Ross et al. 2023). One batch job runs the truth on 256², 64² without closure
   and 64² with the Reynolds backscatter at the published C_R = 7 (10 years, 5 seeds each), then
   computes the KE spectrum, the resolved KE and APE transfers and the subgrid energy transfer (closure,
   and the exact subgrid forcing of the coarse-grained truth with a Gaussian and a sharp filter), and
   draws them in the layout of their Fig. 8:
   ```bash
   sbatch scripts/perezhogin_validation.sh            # from the QGTestbed root folder
   ```
   Settings (latitude, grids, C_R, years, seeds) are at the top of the script. Runs, logs,
   `spectra.jld2` and `perezhogin_validation.png` go to `perezhogin-validation/`. The two analysis
   steps can be repeated on their own:
   ```bash
   julia --project=. scripts/perezhogin_validation_spectra.jl perezhogin-validation 64
   julia --project=plotting scripts/plot_perezhogin_validation.jl perezhogin-validation
   ```
   The paper's JansenHeld run is not included. Spectra are normalized so that their sum times dk is the
   domain mean, with k the angular wavenumber; if pyqg normalizes differently, absolute values can
   differ from the paper by a constant factor while shapes and ratios stay comparable.

9. **All closures in the same layout** (`scripts/closure_comparison.sh`): truth and, on 64², no
   closure, Reynolds averaging (NEMO default with GM, and `nemo_eorca025`: our eORCA025 namelist a3ue,
   no GM, bilaplacian rn_Uv = 0.0838), Mak 2023, Jansen 2019 and Perezhogin 2023 (C_R = 7), as
   published, at 40°N by default (20 years, spin-up 8, 10 seeds, sweep numerics: no dealiasing,
   `QG_FILTER_ALL=true`). Panel (d) shows each closure's full energy transfer (GM, viscosity and
   backscatter; the exponential filter is not included) against the exact subgrid transfer of the
   truth. For Mak and Jansen the eddy-energy field is stored with each PV snapshot
   (`closure_snapshots` in the run file) so their tendency can be recomputed offline.
   ```bash
   sbatch scripts/closure_comparison.sh                               # 40°N
   LAT=30 sbatch scripts/closure_comparison.sh                        # another latitude
   julia --project=. scripts/closure_comparison_spectra.jl closure-comparison/40N_nx64 64
   julia --project=plotting scripts/plot_closure_comparison.jl closure-comparison/40N_nx64 [none,nemo,reynolds]
   ```
   Output in `closure-comparison/<lat>N_nx64/`: `closure_comparison.png` (all closures),
   `closure_comparison_none-nemo-reynolds.png`, the CSV curves and a table of KE/truth and the
   integrated transfers in `logs/spectra.log`.

---

## 7b. Tuning exercise (separate from the as-published comparison)

The comparison in Sections 4–7 uses every closure with its published settings. This exercise
asks the follow-up question: **with its best parameters, how close does each closure get to the
truth, and do the best parameters stay the same across Δx/Ld?** A closure whose optimum does not
move with Δx/Ld is scale aware in practice; one that needs different values at different
latitudes or grids is not.

Tuning runs go to `output/tune/<name>/`, with the parameter values in the file name and stored in
the file (`metadata`), so they never overwrite the as-published sweep in `output/`. Every candidate
parameter set is a short coarse run (10 years, first 3 excluded) scored against the truth with the
same metrics as `scale_awareness.jl`. All candidates run in parallel.

### Prerequisites

- Truth runs (Section 4) for every latitude you tune or validate on.
- Optional: the as-published coarse runs (Section 5) of the same cases. `tune_evaluate.jl` then
  shows their cost as the reference line.

### Tunable parameters

Names as typed on the command line, in convenient units (converted to SI by `scripts/tune_params.jl`):

| closure | parameters | published values |
|---|---|---|
| `nemo` | `Uv` (m/s, bilaplacian), `Ue` (m/s, GM), `Le_km` | Uv = 0.1, Ue = 0.02, Le = Δx |
| `geometric`, `geometric_nosplit` | `alpha`, `lambda_days`, `etaE` (m²/s), `E0` (m³/s²), `L_split_km`, `Uv` | 0.06, 80, 500, 5, 100, 0.5 |
| `meke` | `c`, `c_smag`, `tau_vis_days`, `c_d`, `U_bg` (m/s) | 0.15, 0.06, 45, 0.003, 0.1 |
| `reynolds` | `C_R`, `FGR` | 7, 2 |

### Step 1 — calibration runs

Calibrate on two gray-zone cases: 40°N on 64² (Δx/Ld = 0.89) and 60°N on 64² (Δx/Ld = 1.20),
with two seeds. Each grid contains the published value, so the published setting is part of the
comparison. One job list for all four closures:

```bash
{
julia scripts/tune_plan.jl --closure reynolds  --param C_R=0,1,3,7,15,30 \
      --cases 40:64,60:64 --seeds 1,2
julia scripts/tune_plan.jl --closure nemo      --param Uv=0.025,0.05,0.1 --param Ue=0,0.005,0.01,0.02,0.04 \
      --cases 40:64,60:64 --seeds 1,2
julia scripts/tune_plan.jl --closure geometric --param alpha=0.02,0.04,0.06,0.1 --param lambda_days=20,80,320 \
      --cases 40:64,60:64 --seeds 1,2
julia scripts/tune_plan.jl --closure meke      --param c=0.05,0.1,0.15,0.3,0.6 \
      --cases 40:64,60:64 --seeds 1,2
} > tune_jobs.txt
```

This gives 152 runs (reynolds 24, nemo 60, geometric 48, meke 20). Options of `tune_plan.jl`:
`--closure`, `--param NAME=v1,v2,…` (repeat for more parameters; all combinations are run),
`--cases lat:grid,…`, `--years` (default 10), `--spinup` (default 3), `--seeds` (default 1),
`--outdir` (default `output/tune/<closure>`).

### Step 2 — run

```bash
sbatch scripts/slurm/tune.sh                     # reads tune_jobs.txt
```

**More seeds.** `scripts/slurm/tune_seeds.sh` writes its own job list and runs it, with everything set
by environment variables. Runs go to the same `output/tune/<closure>/` folder, so existing seeds are
reused; e.g. add seeds 2 and 3 to a seed-1 experiment:
```bash
CLOSURE=reynolds PARAMS="C_R=0,1,3,7,15,30" CASES=40:64,60:64 SEEDS=2,3 sbatch scripts/slurm/tune_seeds.sh
CLOSURE=geometric PARAMS="alpha=0.02,0.04,0.06,0.1 lambda_days=20,80,320" SEEDS=1,2,3 sbatch scripts/slurm/tune_seeds.sh
```
`tune_evaluate.jl` averages the cost over all seeds it finds. Defaults: `CASES=40:64,60:64`,
`SEEDS=1,2,3`, `YEARS=10`, `SPINUP=3`.

GNU parallel, 100 runs at a time, single-threaded. Resubmit if the job reaches its time limit:
successful runs are skipped, failed ones are rerun. Per-run logs are in `logs/tune/`.

### Step 3 — score and pick the best set

```bash
for cl in reynolds nemo geometric meke; do
    julia --project=plotting scripts/tune_evaluate.jl --closure $cl
done
```

For each closure this prints the five best parameter sets (mean cost over cases and seeds) and the
cost of the published setting, and writes into `output/tune/<closure>/`:
- `tune_results.csv`: every parameter set ranked, with cost and the three metrics per case and as mean;
- `tune_cost.png`: cost against the parameter (one parameter, one line per case, dashed = published)
  or a cost map (two parameters). A large gold star marks the best set (mean over cases), small stars
  the optimum of each case; a box beside the plot lists the best values, their cost per case next to the
  published setting, and each case's own optimum.

Cost = w₁·spectral error + w₂·|resolved KE error| + w₃·|PV-flux error|, default weights 1,1,1.
Change the balance with `--weights`, e.g. `--weights 1,1,0` to ignore the PV flux. Other options:
`--dir`, `--truthdir` (default `output`), `--spinup` (default 3), `--truth-spinup` (default 5).

Look at the per-case rows too: if the best value at 40°N differs clearly from the best at 60°N,
the closure is already not scale aware between these two cases.

### Step 4 — validation across the sweep

Run each closure's best set (replace the values below with your results from Step 3) on all
other cases, into a separate folder:

```bash
VAL="20:32,20:64,20:128,30:32,30:64,30:128,40:32,40:128,50:32,50:64,50:128,60:32,60:128,70:32,70:64,70:128"
{
julia scripts/tune_plan.jl --closure reynolds  --param C_R=<best>                          --cases $VAL --outdir output/tune/reynolds_validation
julia scripts/tune_plan.jl --closure nemo      --param Uv=<best> --param Ue=<best>        --cases $VAL --outdir output/tune/nemo_validation
julia scripts/tune_plan.jl --closure geometric --param alpha=<best> --param lambda_days=<best> --cases $VAL --outdir output/tune/geometric_validation
julia scripts/tune_plan.jl --closure meke      --param c=<best>                            --cases $VAL --outdir output/tune/meke_validation
} > validation_jobs.txt
JOBS=validation_jobs.txt sbatch scripts/slurm/tune.sh
```

Score them the same way:

```bash
for cl in reynolds nemo geometric meke; do
    julia --project=plotting scripts/tune_evaluate.jl --closure $cl --dir output/tune/${cl}_validation
done
```

The per-case rows of `tune_results.csv` show how the tuned closure performs at every Δx/Ld,
next to the cost of the published setting at the same case.

### Interpretation and caveats

- Compare per case: tuned vs published cost, and how the tuned cost changes with Δx/Ld. A tuned set
  that stays good across all cases indicates a scale-aware closure.
- Tuning runs are shorter than the sweep runs (10 vs 20 years), so their time means are noisier.
  Use at least two seeds and don't over-interpret small cost differences.
- The cost depends on the weights. Report which weights were used, and check that the ranking of the
  closures doesn't hinge on them (e.g. repeat with `--weights 1,1,0` and `--weights 0,0,1`).
- `nemo`, `geometric` and `geometric_nosplit` include NEMO's bilaplacian viscosity (`Uv`); tuning `Uv`
  changes the dissipation package, not the eddy closure itself.

## 7c. Repeating the sweep with the tuned parameters

`scripts/slurm/tuned_sweep.sh` reruns the whole coarse sweep with the best values from Section 7b.
Edit the settings block at the top of the script, one line per closure:

```bash
NEMO_PARAMS="Uv=0.1 Ue=0.02"
GEOMETRIC_PARAMS="alpha=0.06 lambda_days=80"
GEOMETRIC_NOSPLIT_PARAMS="alpha=0.06 lambda_days=80"
MEKE_PARAMS="c=0.15"
REYNOLDS_PARAMS="C_R=3"

CLOSURES="none nemo geometric geometric_nosplit meke reynolds"
LATS="20 30 40 50 60 70";  GRIDS="32 64 128";  SEEDS="1 2 3";  YEARS=20;  SPINUP=5
OUTDIR=output/tuned
```

Parameters not listed keep their published values; an empty string runs the closure as published;
remove a closure from `CLOSURES` to skip it (`none` is included as the reference). Then:

```bash
sbatch scripts/slurm/tuned_sweep.sh
```

The script checks the parameter names before starting, records the settings in
`output/tuned/tuned_settings.txt` (they are also stored in every output file as `metadata`), and runs
the sweep with GNU parallel (resubmit to finish an interrupted sweep). The files have the same names
as in `output/`, which stays untouched. Compare with the truth runs in `output/`:

```bash
julia --project=plotting scripts/scale_awareness.jl -d output/tuned -t output
```

This writes the scale-awareness table and figures into `output/tuned/`; put them next to the
as-published ones in `output/` to see what tuning changed at each Δx/Ld.

## 7d. Numerics: dealiasing, the pyqg validation and the stability test

**Dealiasing is off by default (`aliased_fraction = 0`).** GeophysicalFlows applies the 2/3 rule to the
nonlinear term by default; pyqg does not and relies on its exponential filter alone. With the 2/3 rule,
the 64² runs had about half of pyqg's eddy energy transfer, so every coarse model, with or without a
closure, was too dissipative in the gray zone. Without it, the testbed reproduces pyqg and
Perezhogin & Glazunov (2023, Fig. 10) at 64²: the KE spectra, the resolved KE and APE transfers and the
subgrid transfer of the Reynolds closure agree with pyqg running the authors' own `Reynolds_stress`
code within a few percent (ensemble scatter). The validation is reproducible with
`scripts/perezhogin_validation.sh` (testbed) and the separate `pyqg-validation` package (pyqg, with
`analyze_pyqg.py --testbed` overlaying both).

The setting is read from the environment variable `QG_ALIASED_FRACTION` (default 0); all batch scripts
export it and print it, and every output file stores `aliased_fraction` and `stepper`.
**Runs made before this change used 1/3** (files without an `aliased_fraction` entry). Do not mix them
with new runs in one analysis; to reproduce them, submit with `QG_ALIASED_FRACTION=0.3333333`.

**Stability test.** `nemo`, `geometric`, `geometric_nosplit` and `meke` bring their own dissipation and
run with plain RK4, without the exponential filter. Without dealiasing, nothing else removes grid-scale
aliasing, so check them before the sweep:

```bash
sbatch scripts/slurm/stability.sh                       # closures as published (plain RK4)
cat output/stability_filter-false/stability_summary.txt
```

All closures on 64² at all six latitudes, 5 years each (36 runs, about an hour); set `GRIDS="32 64 128"` (and `LATS`) when submitting for more.
`scripts/stability_check.jl` flags runs that stopped on NaN (`UNSTABLE`), pile up energy above 0.65 k_max
(`NOISY`: grid fraction > 0.02 and > 3× the filtered no-closure run) or keep growing (`GROWING`).
If all runs are `ok`, keep `QG_FILTER_ALL=false`. Otherwise repeat with the filter for every closure:

```bash
QG_FILTER_ALL=true sbatch scripts/slurm/stability.sh
cat output/stability_filter-true/stability_summary.txt
```

and, if that is clean, submit every job of the sweep with `QG_FILTER_ALL=true` (as in pyqg, where all
runs have the filter; report it as a numerical choice, the closures' own dissipation is unchanged).

### Rerunning everything with the new numerics

```bash
# 0. keep the old (dealiased) results apart
mv output output_dealiased && mv logs logs_dealiased
rm -f *_joblog.txt                     # old GNU parallel job logs would mark runs as done
# 1. update the code, then on the login node
julia setup.jl && julia setup.jl plotting
# 2. stability test (decides QG_FILTER_ALL), see above
# 3. truth runs, coarse sweep, analysis (Sections 4, 5, 7), tuning (7b, 7c) as before,
#    with QG_FILTER_ALL=true in front of every sbatch if step 2 requires it
```

---

## 8. How each closure is implemented (2-layer translation)

All closures add a PV tendency per layer through the `calcFq` hook of `MultiLayerQG`.
Prognostic eddy energy is stepped once per time step.

| key | time stepper / base dissipation | eddy part |
|---|---|---|
| `none` | FilteredRK4 (exponential spectral filter, as in pyqg) | – |
| `nemo` | RK4, bilaplacian ahm = U_v Δx³/12 | GM, κ = U_e L_e/2 |
| `geometric(_nosplit)` | RK4, bilaplacian as `nemo` | GM, κ = α Ê/∫(M²/N)dz with prognostic Ê |
| `meke` | RK4, biharmonic Smagorinsky | GM + negative Laplacian viscosity ν₂ = −c√(2e)L_mix (not scaled by R, as in the authors' MOM6 code), prognostic MEKE |
| `reynolds` | FilteredRK4 (as in the paper's pyqg runs) | −C_R ∇·σ^KEB on PV, both layers |

- GM in QG = down-gradient diffusion of the stretching part of PV (interface height):
  ∂q₁ = −F₁∇·(κ∇Δψ), ∂q₂ = +F₂∇·(κ∇Δψ), Δψ = ψ₁ − ψ₂ including the imposed shear.
- GEOMETRIC in 2 layers: κ = α Ê √(g′/H)/(f₀|Δu|), source κ f₀²|Δu|²/g′ (from NΔz = √(g′H)).
- Momentum viscosities act as the curl of a momentum tendency, so the local dissipation that feeds
  MEKE integrates exactly to the resolved energy loss.
- The energy couplings (GM, Laplacian, biharmonic) were verified numerically: the resolved energy
  loss equals the energy given to the eddy-energy budget to round-off.

### Settings: published values, no retuning

| closure | values used | source |
|---|---|---|
| `nemo` | U_v = 0.1 m/s, U_e = 0.02 m/s, L_e = Δx | NEMO 4 grid-scaled coefficients; **set to your namelist** (namdyn_ldf, namtra_eiv) |
| `geometric` | α = 0.06, λ⁻¹ = 80 d, η_E = 500 m²/s, Ê₀ = 5 m³/s², 10 ≤ κ ≤ 3·10⁴ m²/s, splitting L = 100 km, M = 2, bilaplacian U_v = 0.5 m/s | Mak et al. (2023): the authors' NEMO 4.0.5 code and namelists (R025/R050); see below |
| `meke` | c = 0.15, c_smag = 0.06, τ_vis = 45 d, c_d = 0.003, U_bg = 0.1 m/s | Jansen et al. (2019): the authors' MOM6 parameter files ([MFJansen/MEKE_GM_BS_parameters](https://github.com/MFJansen/MEKE_GM_BS_parameters)); see below |
| `reynolds` | C_R = 7, filter-to-grid ratio 2, Gaussian filter | Perezhogin & Glazunov (2023), QG experiments |

Choices where a paper leaves something open or the 2-layer translation needs one:
- Jansen et al.: 1/H in the bottom-drag term for consistent units (as in MOM6); backscatter with the
  same ν₂ in both layers (vertical structure not specified in the paper).
- Mak et al.: the bilaplacian formula is used; the values in their table are about 5× larger.
- Perezhogin & Glazunov: `dynamic_CR` exists in the code but is not part of the paper and is off.

### Jansen et al. (2019): paper text vs the runs that produced the paper

The MEKE closure follows the authors' published MOM6 parameter files and MOM6 code, which differ
from the paper's text in four places. The parameter files are what was run, so they are the default;
the paper-text variant remains available through closure options:

| | paper text | authors' runs (default here) | option |
|---|---|---|---|
| drag coefficient | c_d = 0.03 | `CDRAG = 0.003` | `c_d` |
| backscatter viscosity | ν₂ = −K_GM (includes R) | `MEKE%Ku = −0.15·√(2e)·L_mix`, added without R (`MOM_hor_visc`) | `bs_resolution_scaled` |
| drag velocity | √(e + u_bot² + U_bg²) | √(2·MEKE + u_bot² + `DRAG_BG_VEL`²), `MEKE_CD_SCALE = 1` | `drag_eke_factor` (2 or 1) |
| biharmonic viscosity | (c_smag·\|D\| + 1/τ)·Δ⁴ | max(Smagorinsky, Δ⁴/`AH_TIME_SCALE`) | `nu4_max` |

Unchanged and confirmed by the files: c = 0.15 (`MEKE_KHCOEFF`), R = 1/(1 + π·Ld/Δx)
(`KH_RES_SCALE_COEF = π`, `KH_RES_FN_POWER = 1`), L_mix = min(Δ, Rhines scale), K_e = c·√(2e)·L_mix,
c_smag = 0.06, τ_vis = 45 d, U_bg = 0.1 m/s, drag divided by the column depth.
Runs with `meke` made before this change used the paper-text variant.

### Mak et al. (2023): checked against the authors' NEMO code

`ldfeke.F90` and `ldftra_split.F90` of the authors' NEMO 4.0.5 code implement the same formulas
as the testbed: κ_gm = α·Ê / ∫(M²/N) dz, source ∫κ·M⁴/N² dz, dissipation −λ(Ê − Ê₀), advection by
the depth-mean flow and the long Rossby speed, Laplacian diffusion of Ê, and splitting by solving
(1 − L²∇²) twice on T and S, with slopes and N² then taken from the large-scale field. The namelists of
their eddy-permitting runs (R010–R050) differ from the earlier defaults here in three values, now
adopted: bilaplacian viscosity `rn_Uv = 0.5` m/s (was 0.1), `rn_eke_min = 5` m³/s² (was 4, as in the
paper's table) and the bounds 10 ≤ κ_gm ≤ 3·10⁴ m²/s (was uncapped). Their taper of κ_gm for
deformation radii below 7 km is inactive at all testbed latitudes and not reproduced. Their R010
run uses a 200 km splitting scale; the testbed keeps the 100 km of R025 and R050.
Runs with `geometric*` made before this change used the earlier values.

### Known limitations of this testbed

- The mean shear is imposed and never changes, so GM cannot flatten the mean interface. GM-type
  closures act only on the resolved eddies; their mean PV flux is diagnosed (`vq_param`), not applied.
- With the published Ê₀ and λ, GEOMETRIC's eddy energy in this configuration settles near
  λÊ₀/(λ − α f₀ΔU/√(g′H)), i.e. it is set mainly by Ê₀ and λ.
- Ld varies only by a factor of about 2.8 across the sweep (realistic ocean: more). To widen it,
  make g′ latitude dependent in `src/cases.jl`.
- QG is invalid near the equator; the sweep starts at 20°N.

---

## 9. Repository layout

```
setup.jl                      create Project.toml, install model (and plotting) environments, precompile
src/cases.jl                  latitude-dependent parameters
src/problem.jl                MultiLayerQG problem + closure hook
src/closures/*.jl             the closures
src/run.jl                    run_simulation (driver and output)
src/diagnostics.jl            energies, spectra, PV fluxes
src/coarsegrain.jl            filtering, subgrid forcing, a priori skill
scripts/smoke_test.jl         installation check
scripts/run_member.jl         one run from the command line
scripts/sweep.jl              writes the full job list
scripts/apriori.jl            offline skill on truth snapshots
scripts/plot_vorticity.jl     vorticity comparison panels
scripts/plot_summary.jl       energy time series and KE spectra
scripts/plot_circulation.jl   mean flow and eddies: SSH contours, transport shading (Jansen et al. Fig. 3)
scripts/scale_awareness.jl    metrics of every run vs truth, plotted against Ld/Δx
scripts/meke_budget.jl        energy budget of the MEKE closure (what sets K_GM)
scripts/stability_check.jl    stability test of the closures without dealiasing (Section 7d)
scripts/slurm/stability.sh    Slurm job: stability test (Section 7d)
scripts/perezhogin_validation.sh         Slurm job: Perezhogin & Glazunov (2023) Fig. 8 validation at 49°N
scripts/perezhogin_validation_run.jl     one run of that validation (truth / none / reynolds)
scripts/perezhogin_validation_spectra.jl its spectra and energy transfers (model environment)
scripts/plot_perezhogin_validation.jl    its figure (plotting environment)
scripts/closure_comparison.sh            Slurm job: all closures in the Fig. 10 layout (default 40°N, 64²)
scripts/closure_comparison_run.jl        one run of that comparison (truth or one closure)
scripts/closure_comparison_spectra.jl    its spectra and energy transfers, incl. each closure's transfer
scripts/plot_closure_comparison.jl       its figure (plotting environment; optional closure subset)
scripts/plotting_utils.jl     shared spectrum helpers for the plotting scripts
scripts/slurm/truth.sh        Slurm job: truth runs
scripts/slurm/truth_and_50N.sh  Slurm job: truth runs + 50°N 64² closure runs
scripts/slurm/coarse.sh       Slurm job: coarse runs
scripts/slurm/tune.sh         Slurm job: tuning runs
scripts/slurm/tune_seeds.sh   Slurm job: tuning runs with several seeds (writes its own job list)
scripts/slurm/tuned_sweep.sh  Slurm job: full sweep with the tuned parameters (settings block at the top)
scripts/tune_plan.jl          tuning: write the job list for a parameter grid
scripts/tune_member.jl        one run with given closure parameters (tuning and tuned sweep)
scripts/tune_evaluate.jl      tuning: score against the truth, rank, plot
scripts/tune_params.jl        tuning: tunable parameter names per closure
```

---

## 10. Troubleshooting

| message | cause | fix |
|---|---|---|
| `Package QGTestbed does not have LinearAlgebra in its dependencies` | `Project.toml` lost its dependencies (environment created without `setup.jl`, or overwritten by an older package version that shipped a bare `Project.toml`) | `rm Project.toml Manifest.toml`, then `julia setup.jl` on the login node |
| `QGTestbed Being precompiled by another process` in many job logs | jobs started before QGTestbed was compiled | run `julia setup.jl` on the login node before submitting (the batch scripts now also compile once before starting the runs) |
| `Package JLD2 not found in current path` when plotting | plotting script run outside the plotting environment | `julia setup.jl plotting`, then `julia --project=plotting …` from this folder |
| `Circular dependency detected` / `Failed to precompile Makie` | Makie cannot be precompiled with Julia 1.11.0 (extension cycle) | the plotting environment uses Plots.jl instead; if an old `plotting/` folder with Makie exists, delete it and run `julia setup.jl plotting` again |
| batch job ends immediately without output | the GNU parallel job log (`coarse_joblog.txt`, `tune_jobs_joblog.txt`) lists all jobs from an earlier attempt; with `--resume` they are all skipped | delete the job log and resubmit (the scripts now use `--resume-failed`, which reruns failed jobs) |
| `Assignment to prob in soft scope is ambiguous` | Julia scoping warning in scripts | harmless |
| runs much slower than expected | FFTs planned for all cores of the node | set `JULIA_CPU_THREADS` equal to `-t` (Section 4) |

---

## Appendix — batch scripts

All scripts are in `scripts/slurm/` and are submitted from the `QGTestbed` folder with `sbatch`.
Before the first use, replace `<your_account>` and `<your_qos>`, and check that `module load julia` loads the
same Julia version that was used for `julia setup.jl`. Each script first loads QGTestbed once (so the runs don't
all compile it at the same time) and stops with a message if the environment is broken. Runs started with `&`
inside a script must be followed by `wait`, otherwise the job ends immediately and Slurm stops the runs.

### Truth runs (Section 4): `scripts/slurm/truth.sh`

```bash
#!/bin/bash
# Truth (eddy-resolving, no closure) runs for all latitudes, one seed, on one node.
# Submit from the QGTestbed folder:  sbatch scripts/slurm/truth.sh
#SBATCH --job-name=qg_truth
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --time=48:00:00
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=truth_%j.out

module load julia                  # the same Julia module used for `julia setup.jl`
# Numerics (README section 7d): no 2/3 dealiasing (validated against pyqg); exponential filter for
# all closures only if QG_FILTER_ALL=true (decide with scripts/slurm/stability.sh). Override when
# submitting, e.g.  QG_FILTER_ALL=true sbatch ...
export QG_ALIASED_FRACTION=${QG_ALIASED_FRACTION:-0}
export QG_FILTER_ALL=${QG_FILTER_ALL:-false}
echo "numerics: aliased_fraction = $QG_ALIASED_FRACTION, exponential filter for all closures = $QG_FILTER_ALL"

# Compile QGTestbed once before starting the runs, so they don't all compile it at the same time
julia --project=. -e 'using QGTestbed' || { echo "QGTestbed does not load: run 'julia setup.jl' on the login node"; exit 1; }

SEED=${SEED:-1}
YEARS=${YEARS:-20}
SPINUP=${SPINUP:-5}
mkdir -p logs

# 512² grids (50–70°N): 16 FFT threads each. JULIA_CPU_THREADS sets the threads
# FourierFlows plans its FFTs for; -t gives Julia the same number.
for lat in 50 60 70; do
    JULIA_CPU_THREADS=16 julia -t 16 --project=. scripts/run_member.jl $lat 0 truth $YEARS $SPINUP $SEED \
        > logs/truth_${lat}N_seed${SEED}.log 2>&1 &
done
# 256² grids (20–40°N): 8 threads each.
for lat in 20 30 40; do
    JULIA_CPU_THREADS=8 julia -t 8 --project=. scripts/run_member.jl $lat 0 truth $YEARS $SPINUP $SEED \
        > logs/truth_${lat}N_seed${SEED}.log 2>&1 &
done
wait
```

### Truth runs plus the 50°N closure runs: `scripts/slurm/truth_and_50N.sh`

```bash
#!/bin/bash
# Truth runs for all latitudes plus the six 64² closure runs at 50°N (for a first
# vorticity comparison before the full coarse sweep), seed 1, on one node.
# Submit from the QGTestbed folder:  sbatch scripts/slurm/truth_and_50N.sh
# Afterwards remove the 50°N 64² seed-1 runs from coarse_jobs.txt (see README, section 5).
#SBATCH --job-name=qg_truth50
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --time=48:00:00
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=truth50_%j.out

module load julia                  # the same Julia module used for `julia setup.jl`
# Numerics (README section 7d): no 2/3 dealiasing (validated against pyqg); exponential filter for
# all closures only if QG_FILTER_ALL=true (decide with scripts/slurm/stability.sh). Override when
# submitting, e.g.  QG_FILTER_ALL=true sbatch ...
export QG_ALIASED_FRACTION=${QG_ALIASED_FRACTION:-0}
export QG_FILTER_ALL=${QG_FILTER_ALL:-false}
echo "numerics: aliased_fraction = $QG_ALIASED_FRACTION, exponential filter for all closures = $QG_FILTER_ALL"

# Compile QGTestbed once before starting the runs, so they don't all compile it at the same time
julia --project=. -e 'using QGTestbed' || { echo "QGTestbed does not load: run 'julia setup.jl' on the login node"; exit 1; }
mkdir -p logs

# 512² truth runs: 16 threads each (48 cores)
for lat in 50 60 70; do
    JULIA_CPU_THREADS=16 julia -t 16 --project=. scripts/run_member.jl $lat 0 truth 20 5 1 > logs/truth_${lat}N.log 2>&1 &
done
# 256² truth runs: 8 threads each (24 cores)
for lat in 20 30 40; do
    JULIA_CPU_THREADS=8 julia -t 8 --project=. scripts/run_member.jl $lat 0 truth 20 5 1 > logs/truth_${lat}N.log 2>&1 &
done
# 64² closure runs at 50°N: 1 thread each (6 cores)
for cl in none nemo geometric geometric_nosplit meke reynolds; do
    JULIA_CPU_THREADS=1 julia -t 1 --project=. scripts/run_member.jl 50 64 $cl 20 5 1 > logs/${cl}_50N.log 2>&1 &
done

wait    # keep the job alive until every run above has finished
```

### Coarse sweep (Section 5): `scripts/slurm/coarse.sh`

```bash
#!/bin/bash
# All coarse runs (latitudes × grids × closures × seeds), single-threaded, 100 at a time.
# Create the job list first, then submit from the QGTestbed folder:
#   julia --project=. scripts/sweep.jl | grep -v truth | sed 's/^/JULIA_CPU_THREADS=1 /' > coarse_jobs.txt
#   sbatch scripts/slurm/coarse.sh
# If the job hits the time limit, submit it again: --resume-failed skips successful runs and reruns failed ones.
#SBATCH --job-name=qg_coarse
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --time=48:00:00
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=coarse_%j.out

module load julia                  # the same Julia module used for `julia setup.jl`
# Numerics (README section 7d): no 2/3 dealiasing (validated against pyqg); exponential filter for
# all closures only if QG_FILTER_ALL=true (decide with scripts/slurm/stability.sh). Override when
# submitting, e.g.  QG_FILTER_ALL=true sbatch ...
export QG_ALIASED_FRACTION=${QG_ALIASED_FRACTION:-0}
export QG_FILTER_ALL=${QG_FILTER_ALL:-false}
echo "numerics: aliased_fraction = $QG_ALIASED_FRACTION, exponential filter for all closures = $QG_FILTER_ALL"
module load parallel 2>/dev/null   # GNU parallel, if not on the default path

# Compile QGTestbed once before starting the runs, so they don't all compile it at the same time
julia --project=. -e 'using QGTestbed' || { echo "QGTestbed does not load: run 'julia setup.jl' on the login node"; exit 1; }
export JULIA_NUM_THREADS=1

mkdir -p logs
parallel -j 100 --joblog coarse_joblog.txt --resume-failed --results logs/ < coarse_jobs.txt
```

### Tuning runs (Section 7b): `scripts/slurm/tune.sh`

```bash
#!/bin/bash
# Tuning runs from a job list written by scripts/tune_plan.jl, single-threaded, 100 at a time.
#   julia scripts/tune_plan.jl --closure reynolds --param C_R=1,3,7,15,30 --cases 40:64,50:64 > tune_jobs.txt
#   sbatch scripts/slurm/tune.sh                 # or: JOBS=other_list.txt sbatch scripts/slurm/tune.sh
# Resubmit to finish interrupted lists: --resume-failed skips successful runs and reruns failed ones.
#SBATCH --job-name=qg_tune
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --time=48:00:00
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=tune_%j.out

module load julia                  # the same Julia module used for `julia setup.jl`
# Numerics (README section 7d): no 2/3 dealiasing (validated against pyqg); exponential filter for
# all closures only if QG_FILTER_ALL=true (decide with scripts/slurm/stability.sh). Override when
# submitting, e.g.  QG_FILTER_ALL=true sbatch ...
export QG_ALIASED_FRACTION=${QG_ALIASED_FRACTION:-0}
export QG_FILTER_ALL=${QG_FILTER_ALL:-false}
echo "numerics: aliased_fraction = $QG_ALIASED_FRACTION, exponential filter for all closures = $QG_FILTER_ALL"
module load parallel 2>/dev/null

# Compile QGTestbed once before starting the runs, so they don't all compile it at the same time
julia --project=. -e 'using QGTestbed' || { echo "QGTestbed does not load: run 'julia setup.jl' on the login node"; exit 1; }
export JULIA_NUM_THREADS=1

JOBS=${JOBS:-tune_jobs.txt}
mkdir -p logs/tune
parallel -j 100 --joblog ${JOBS%.txt}_joblog.txt --resume-failed --results logs/tune/ < "$JOBS"
```

### Tuning runs with several seeds (Section 7b): `scripts/slurm/tune_seeds.sh`

```bash
#!/bin/bash
# Tuning runs with several seeds: writes its own job list with scripts/tune_plan.jl and runs it.
# All settings come from environment variables, so no job list has to be prepared beforehand.
#
#   CLOSURE=reynolds PARAMS="C_R=0,1,3,7,15,30" CASES=40:64,60:64 SEEDS=2,3 \
#       sbatch scripts/slurm/tune_seeds.sh
#
#   CLOSURE=geometric PARAMS="alpha=0.02,0.04,0.06,0.1 lambda_days=20,80,320" \
#       CASES=40:64,60:64 SEEDS=1,2,3 sbatch scripts/slurm/tune_seeds.sh
#
#   CLOSURE   closure to tune (nemo | geometric | geometric_nosplit | meke | reynolds)   [required]
#   PARAMS    space-separated NAME=v1,v2,… (one per parameter, all combinations run)    [required]
#   CASES     latitude:grid pairs                                         (default 40:64,60:64)
#   SEEDS     seeds to run                                                (default 1,2,3)
#   YEARS     run length in years                                         (default 10)
#   SPINUP    spin-up excluded from the scores                            (default 3)
#
# Runs go to output/tune/<closure>/ like the single-seed tuning, so seeds already run there
# (e.g. seed 1 from scripts/slurm/tune.sh) are reused: set SEEDS=2,3 to add only new seeds.
# tune_evaluate.jl averages the cost over all seeds it finds.
# Resubmit with the same settings to finish an interrupted list (--resume-failed).
#SBATCH --job-name=qg_tune_seeds
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --time=48:00:00
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=tune_seeds_%j.out

module load julia                  # the same Julia module used for `julia setup.jl`
# Numerics (README section 7d): no 2/3 dealiasing (validated against pyqg); exponential filter for
# all closures only if QG_FILTER_ALL=true (decide with scripts/slurm/stability.sh). Override when
# submitting, e.g.  QG_FILTER_ALL=true sbatch ...
export QG_ALIASED_FRACTION=${QG_ALIASED_FRACTION:-0}
export QG_FILTER_ALL=${QG_FILTER_ALL:-false}
echo "numerics: aliased_fraction = $QG_ALIASED_FRACTION, exponential filter for all closures = $QG_FILTER_ALL"
module load parallel 2>/dev/null

: "${CLOSURE:?set CLOSURE, e.g. CLOSURE=reynolds}"
: "${PARAMS:?set PARAMS, e.g. PARAMS=\"C_R=0,1,3,7,15,30\"}"
CASES=${CASES:-40:64,60:64}
SEEDS=${SEEDS:-1,2,3}
YEARS=${YEARS:-10}
SPINUP=${SPINUP:-3}

# Compile QGTestbed once before starting the runs, so they don't all compile it at the same time
julia --project=. -e 'using QGTestbed' || { echo "QGTestbed does not load: run 'julia setup.jl' on the login node"; exit 1; }
export JULIA_NUM_THREADS=1

# Job list for these settings (one file per closure and seed set)
JOBS=tune_${CLOSURE}_seeds_${SEEDS//,/-}.txt
PARAM_ARGS=()
for p in $PARAMS; do PARAM_ARGS+=(--param "$p"); done
julia scripts/tune_plan.jl --closure "$CLOSURE" "${PARAM_ARGS[@]}" --cases "$CASES" \
      --seeds "$SEEDS" --years "$YEARS" --spinup "$SPINUP" > "$JOBS" || { echo "tune_plan.jl failed"; exit 1; }
echo "$(wc -l < "$JOBS") runs in $JOBS"

mkdir -p logs/tune
parallel -j 100 --joblog ${JOBS%.txt}_joblog.txt --resume-failed --results logs/tune/ < "$JOBS"
```

### Sweep with tuned parameters (Section 7c): `scripts/slurm/tuned_sweep.sh`

```bash
#!/bin/bash
# Repeat the coarse sweep with tuned closure parameters.
# Edit the settings block below (best values from tune_evaluate.jl), then submit from the
# QGTestbed folder:   sbatch scripts/slurm/tuned_sweep.sh
#
# Runs go to $OUTDIR (default output/tuned/) with the same file names as the as-published
# sweep in output/, which stays untouched. The parameter values are stored in every file
# (`metadata`) and in $OUTDIR/tuned_settings.txt. Compare with the truth runs in output/:
#   julia --project=plotting scripts/scale_awareness.jl -d output/tuned -t output
# Resubmit with the same settings to finish an interrupted sweep (--resume-failed).
#SBATCH --job-name=qg_tuned_sweep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --time=48:00:00
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=tuned_sweep_%j.out

# ============================ settings ==========================================
# Best parameter values per closure, space-separated NAME=value (names as in tune_plan.jl:
#   nemo: Uv Ue Le_km | geometric(_nosplit): alpha lambda_days etaE E0 L_split_km Uv
#   meke: c c_smag tau_vis_days c_d U_bg | reynolds: C_R FGR).
# Parameters not listed keep their published values. An empty string runs the closure
# with its published values; remove a closure from CLOSURES to skip it.
NEMO_PARAMS="Uv=0.1 Ue=0.02"
GEOMETRIC_PARAMS="alpha=0.06 lambda_days=80"
GEOMETRIC_NOSPLIT_PARAMS="alpha=0.06 lambda_days=80"
MEKE_PARAMS="c=0.15"
REYNOLDS_PARAMS="C_R=3"

CLOSURES="none nemo geometric geometric_nosplit meke reynolds"   # `none` = reference without closure
LATS="20 30 40 50 60 70"
GRIDS="32 64 128"
SEEDS="1 2 3"
YEARS=20
SPINUP=5
OUTDIR=output/tuned
NPAR=100                                                          # runs at a time
# ================================================================================

module load julia                  # the same Julia module used for `julia setup.jl`
# Numerics (README section 7d): no 2/3 dealiasing (validated against pyqg); exponential filter for
# all closures only if QG_FILTER_ALL=true (decide with scripts/slurm/stability.sh). Override when
# submitting, e.g.  QG_FILTER_ALL=true sbatch ...
export QG_ALIASED_FRACTION=${QG_ALIASED_FRACTION:-0}
export QG_FILTER_ALL=${QG_FILTER_ALL:-false}
echo "numerics: aliased_fraction = $QG_ALIASED_FRACTION, exponential filter for all closures = $QG_FILTER_ALL"
module load parallel 2>/dev/null

params_for() {
    case "$1" in
        nemo)              echo "$NEMO_PARAMS" ;;
        geometric)         echo "$GEOMETRIC_PARAMS" ;;
        geometric_nosplit) echo "$GEOMETRIC_NOSPLIT_PARAMS" ;;
        meke)              echo "$MEKE_PARAMS" ;;
        reynolds)          echo "$REYNOLDS_PARAMS" ;;
        *)                 echo "" ;;
    esac
}

# Compile QGTestbed once before starting the runs, so they don't all compile it at the same time
julia --project=. -e 'using QGTestbed' || { echo "QGTestbed does not load: run 'julia setup.jl' on the login node"; exit 1; }
export JULIA_NUM_THREADS=1

# Check the parameter names before starting hundreds of runs (tune_plan.jl validates them)
for cl in $CLOSURES; do
    P=$(params_for $cl)
    [ -z "$P" ] && continue
    ARGS=(); for p in $P; do ARGS+=(--param "$p"); done
    julia scripts/tune_plan.jl --closure $cl "${ARGS[@]}" > /dev/null 2>&1 \
        || { echo "invalid parameters for $cl: '$P' (see scripts/tune_params.jl)"; exit 1; }
done

# Record the settings and write the job list
mkdir -p "$OUTDIR" logs/tuned
{
    echo "tuned sweep, $(date)"
    for cl in $CLOSURES; do
        P=$(params_for $cl); echo "$cl: ${P:-published values}"
    done
    echo "lats: $LATS | grids: $GRIDS | seeds: $SEEDS | years: $YEARS | spin-up: $SPINUP"
} > "$OUTDIR/tuned_settings.txt"
cat "$OUTDIR/tuned_settings.txt"

JOBS=tuned_sweep_jobs.txt
> "$JOBS"
for lat in $LATS; do for nx in $GRIDS; do for seed in $SEEDS; do for cl in $CLOSURES; do
    echo "JULIA_CPU_THREADS=1 julia -t 1 --project=. scripts/tune_member.jl $lat $nx $cl $YEARS $SPINUP $seed $OUTDIR --no-tag $(params_for $cl)" >> "$JOBS"
done; done; done; done
echo "$(wc -l < "$JOBS") runs in $JOBS"

parallel -j $NPAR --joblog tuned_sweep_joblog.txt --resume-failed --results logs/tuned/ < "$JOBS"
```

### Stability test (Section 7d): `scripts/slurm/stability.sh`

```bash
#!/bin/bash
# Stability test of the closures without 2/3 dealiasing (README section 7d).
# nemo, geometric, geometric_nosplit and meke run with plain RK4 and their own dissipation; with
# dealiasing switched off nothing else removes grid-scale aliasing. This job runs every closure
# (published settings) for a few years on all grids at four latitudes, then
# scripts/stability_check.jl reports runs that blew up (NaN) or pile up energy at the grid scale.
# Submit from the QGTestbed folder:
#   sbatch scripts/slurm/stability.sh                       # closures as published (plain RK4)
#   QG_FILTER_ALL=true sbatch scripts/slurm/stability.sh    # every closure with the exponential filter
# Results go to output/stability_filter-<QG_FILTER_ALL>/, a table to stability.csv in that folder.
#SBATCH --job-name=qg_stability
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --time=06:00:00
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=stability_%j.out

# ============================ settings ==========================================
LATS="20 40 60 70"
GRIDS="32 64 128"
CLOSURES="none nemo geometric geometric_nosplit meke reynolds"   # none = filtered reference
YEARS=5
SPINUP=1
SEED=1
NPAR=100
# ================================================================================

module load julia                  # the same Julia module used for `julia setup.jl`
export QG_ALIASED_FRACTION=${QG_ALIASED_FRACTION:-0}
export QG_FILTER_ALL=${QG_FILTER_ALL:-false}
echo "numerics: aliased_fraction = $QG_ALIASED_FRACTION, exponential filter for all closures = $QG_FILTER_ALL"
module load parallel 2>/dev/null
set -o pipefail

julia --project=. -e 'using QGTestbed' || { echo "QGTestbed does not load: run 'julia setup.jl' on the login node"; exit 1; }
export JULIA_NUM_THREADS=1

OUTDIR=output/stability_filter-${QG_FILTER_ALL}
mkdir -p "$OUTDIR" logs/stability
JOBS=$OUTDIR/stability_jobs.txt
> "$JOBS"
for lat in $LATS; do for nx in $GRIDS; do for cl in $CLOSURES; do
    echo "JULIA_CPU_THREADS=1 julia -t 1 --project=. scripts/tune_member.jl $lat $nx $cl $YEARS $SPINUP $SEED $OUTDIR --no-tag" >> "$JOBS"
done; done; done
echo "$(wc -l < "$JOBS") runs"

parallel -j $NPAR --joblog "$OUTDIR/stability_joblog.txt" --resume-failed --results logs/stability/ < "$JOBS"

julia --project=plotting scripts/stability_check.jl -d "$OUTDIR" -s $SPINUP 2>&1 | tee "$OUTDIR/stability_summary.txt"
```

### Perezhogin & Glazunov (2023) validation (Section 7, item 8): `scripts/perezhogin_validation.sh`

```bash
#!/bin/bash
# Validation against Perezhogin & Glazunov (2023), Fig. 10 (Fig. 8 in the preprint): the pyqg "eddy" configuration
# (49°N member of the testbed), truth on 256², 64² without closure and 64² with the
# Reynolds-stress backscatter (C_R = 7), then spectra, energy transfers and the figure.
# Submit from the QGTestbed root folder:
#   sbatch scripts/perezhogin_validation.sh
# Everything (runs, logs, spectra.jld2, perezhogin_validation.png) goes to $OUTDIR, a subfolder of
# perezhogin-validation/ named after the dealiasing setting, so runs with other settings stay apart.
# Needs the model environment (julia setup.jl) and the plotting environment (julia setup.jl plotting).
# Rerunning skips runs whose output file already exists; delete them to start over.
#SBATCH --job-name=qg_perezhogin_val
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --time=24:00:00
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=perezhogin_validation_%j.out

# ============================ settings ==========================================
ALIASED_FRACTION=0  # 0 = no 2/3 dealiasing, exponential filter only, as in pyqg; 0.3333333 = testbed default
OUTDIR=perezhogin-validation/aliased_fraction_${ALIASED_FRACTION}
LAT=49            # 49°N = the pyqg eddy configuration (Ld = 15 km, beta = 1.5e-11, dU = 0.025 m/s)
NX_TRUTH=256      # as in Ross et al. (2023) / Perezhogin & Glazunov (2023)
NX_COARSE=64
C_R=7             # published value for this setup
YEARS=20          # Perezhogin & Glazunov (2023): averages over years 5-20
SPINUP=5
SNAPSHOT_DAYS=30  # PV snapshots for the spectra, after the spin-up (181 per run)
SEEDS="1 2 3 4 5 6 7 8 9 10"   # ensemble of 10 members, as in the paper
TRUTH_THREADS=8   # 10 truth runs x 8 + 20 coarse runs x 1 = 100 cores
# ================================================================================

module load julia                  # the same Julia module used for `julia setup.jl`
set -o pipefail                    # a failing julia step is not hidden by `tee`

julia --project=. -e 'using QGTestbed' || { echo "QGTestbed does not load: run 'julia setup.jl' on the login node"; exit 1; }
julia --project=plotting -e 'using JLD2, Plots' || { echo "plotting environment missing: run 'julia setup.jl plotting'"; exit 1; }

mkdir -p "$OUTDIR/logs"
{
    echo "Perezhogin & Glazunov (2023) validation, $(date)"
    echo "lat $LAT | truth ${NX_TRUTH}² | coarse ${NX_COARSE}² | C_R $C_R | years $YEARS | spin-up $SPINUP | snapshots every $SNAPSHOT_DAYS d | seeds $SEEDS | aliased_fraction $ALIASED_FRACTION"
} > "$OUTDIR/settings.txt"
cat "$OUTDIR/settings.txt"

run() {   # member nx threads seed
    local member=$1 nx=$2 threads=$3 seed=$4
    local file="$OUTDIR/${LAT}N_nx${nx}_${member}_seed${seed}.jld2"
    if [ -f "$file" ]; then echo "exists, skipped: $file"; return; fi
    JULIA_CPU_THREADS=$threads julia -t $threads --project=. scripts/perezhogin_validation_run.jl \
        $member $seed "$OUTDIR" $YEARS $SPINUP $C_R $LAT $NX_TRUTH $NX_COARSE $SNAPSHOT_DAYS $ALIASED_FRACTION \
        > "$OUTDIR/logs/${member}_nx${nx}_seed${seed}.log" 2>&1
}

# 1) simulations, all in parallel
for seed in $SEEDS; do
    run truth    $NX_TRUTH  $TRUTH_THREADS $seed &
    run none     $NX_COARSE 1              $seed &
    run reynolds $NX_COARSE 1              $seed &
done
wait

missing=0
for seed in $SEEDS; do for m in "truth $NX_TRUTH" "none $NX_COARSE" "reynolds $NX_COARSE"; do
    set -- $m
    f="$OUTDIR/${LAT}N_nx${2}_${1}_seed${seed}.jld2"
    [ -f "$f" ] || { echo "missing $f (see $OUTDIR/logs/${1}_nx${2}_seed${seed}.log)"; missing=1; }
done; done
[ $missing -eq 0 ] || echo "some runs failed; the analysis uses the runs that exist"

# 2) spectra and energy transfers (model environment: needs the closure and the coarse-graining)
JULIA_CPU_THREADS=8 julia -t 8 --project=. scripts/perezhogin_validation_spectra.jl "$OUTDIR" $NX_COARSE \
    2>&1 | tee "$OUTDIR/logs/spectra.log" || { echo "spectra failed"; exit 1; }

# 3) figure (plotting environment)
julia --project=plotting scripts/plot_perezhogin_validation.jl "$OUTDIR" 2>&1 | tee "$OUTDIR/logs/plot.log"
echo "done: $OUTDIR/perezhogin_validation.png"
```

---

## References

- Perezhogin, P. & Glazunov, A. (2023). Subgrid parameterizations of ocean mesoscale eddies based on Germano decomposition. JAMES. arXiv:2304.06789
- Mak, J., Marshall, D. P., Madec, G. & Maddison, J. R. (2023). Scale-awareness in an eddy energy constrained mesoscale eddy parameterization. JAMES. arXiv:2306.08988
- Jansen, M. F., Adcroft, A., Khani, S. & Kong, H. (2019). Toward an energetically consistent, resolution aware parameterization of ocean mesoscale eddies. JAMES.
- Ross, A., Li, Z., Perezhogin, P., Fernandez-Granda, C. & Zanna, L. (2023). Benchmarking of machine learning ocean subgrid parameterizations in an idealized model. JAMES.
- Constantinou, N. C. et al. (2021). GeophysicalFlows.jl: Solvers for geophysical fluid dynamics problems in periodic domains on CPUs & GPUs. JOSS 6(60), 3053.
