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
- `Project.toml` is shipped without dependencies; `setup.jl` adds them. Don't run a partial
  `Pkg.add` instead (see Troubleshooting).

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
again: `--resume` with the job log skips completed runs. Per-run logs are in `logs/`.

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
| `snapshots`, `snaptimes` | full PV fields every 30 days after spin-up (truth runs only) |
| `q_final` | final PV field |

Sizes: a few MB for coarse runs, about 250 MB for 512² truth runs.

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
   ```
   Writes `output/vorticity_50N_nx64.png`. Columns: truth (coarse-grained to the same grid)
   and each closure; rows: times; one colour scale per row from the truth.
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
   truth at the same latitude and plot them against Δx/Ld:
   ```bash
   julia --project=plotting scripts/scale_awareness.jl              # all runs in output/
   julia --project=plotting scripts/scale_awareness.jl -s 5 -d output
   ```
   - spectral error: RMS of log₁₀(E_model/E_truth) for 0 < k ≤ k_c/2 (k_c = dealiasing cutoff);
     0 is perfect, 1 means a factor 10 off on average, capped at 6 (eddies gone);
   - resolved KE: model KE / truth KE at k ≤ k_c, minus 1;
   - PV flux: (resolved + parametrized upper-layer ⟨vq⟩) / truth ⟨vq⟩, minus 1.

   Writes `output/scale_awareness.csv` (one row per latitude, grid and closure, mean and range over
   seeds) and one figure per metric, each with one row per grid (32², 64², 128²) and the same Δx/Ld
   axis in every row (colours = closures): `output/scale_awareness_spectral_error.png`,
   `output/scale_awareness_resolved_KE.png`, `output/scale_awareness_PV_flux.png`. A scale-aware
   closure stays near 0 across the whole Δx/Ld range; a non-scale-aware one drifts with Δx/Ld.

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
| `geometric`, `geometric_nosplit` | `alpha`, `lambda_days`, `etaE` (m²/s), `E0` (m³/s²), `L_split_km`, `Uv` | 0.06, 80, 500, 4, 100, 0.1 |
| `meke` | `c`, `c_smag`, `tau_vis_days`, `c_d`, `U_bg` (m/s) | 0.15, 0.06, 45, 0.03, 0.1 |
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

GNU parallel, 100 runs at a time, single-threaded. Resubmit if the job reaches its time limit:
completed runs are skipped. Per-run logs are in `logs/tune/`.

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
  or a cost map (two parameters, star = best).

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

## 8. How each closure is implemented (2-layer translation)

All closures add a PV tendency per layer through the `calcFq` hook of `MultiLayerQG`.
Prognostic eddy energy is stepped once per time step.

| key | time stepper / base dissipation | eddy part |
|---|---|---|
| `none` | FilteredRK4 (exponential spectral filter, as in pyqg) | – |
| `nemo` | RK4, bilaplacian ahm = U_v Δx³/12 | GM, κ = U_e L_e/2 |
| `geometric(_nosplit)` | RK4, bilaplacian as `nemo` | GM, κ = α Ê/∫(M²/N)dz with prognostic Ê |
| `meke` | RK4, biharmonic Smagorinsky | GM + negative Laplacian viscosity ν₂ = −K_GM, prognostic MEKE |
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
| `geometric` | α = 0.06, λ⁻¹ = 80 d, η_E = 500 m²/s, Ê₀ = 4 m³/s², splitting L = 100 km, M = 2, no κ cap | Mak et al. (2023), Table 1 and Sec. 3 |
| `meke` | c = 0.15, c_smag = 0.06, τ_vis = 45 d, c_d = 0.03, U_bg = 0.1 m/s | Jansen et al. (2019) |
| `reynolds` | C_R = 7, filter-to-grid ratio 2, Gaussian filter | Perezhogin & Glazunov (2023), QG experiments |

Choices where a paper leaves something open or the 2-layer translation needs one:
- Jansen et al.: 1/H in the bottom-drag term for consistent units (as in MOM6); backscatter with the
  same ν₂ in both layers (vertical structure not specified in the paper).
- Mak et al.: the bilaplacian formula is used; the values in their table are about 5× larger.
- Perezhogin & Glazunov: `dynamic_CR` exists in the code but is not part of the paper and is off.

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
setup.jl                      install model (and plotting) environments
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
scripts/scale_awareness.jl    metrics of every run vs truth, plotted against Δx/Ld
scripts/plotting_utils.jl     shared spectrum helpers for the plotting scripts
scripts/slurm/truth.sh        Slurm job: truth runs
scripts/slurm/truth_and_50N.sh  Slurm job: truth runs + 50°N 64² closure runs
scripts/slurm/coarse.sh       Slurm job: coarse runs
scripts/slurm/tune.sh         Slurm job: tuning runs
scripts/tune_plan.jl          tuning: write the job list for a parameter grid
scripts/tune_member.jl        tuning: one run with given parameters
scripts/tune_evaluate.jl      tuning: score against the truth, rank, plot
scripts/tune_params.jl        tuning: tunable parameter names per closure
```

---

## 10. Troubleshooting

| message | cause | fix |
|---|---|---|
| `Package QGTestbed does not have LinearAlgebra in its dependencies` | environment created without `setup.jl` | `rm Manifest.toml`, restore `Project.toml` to name/uuid/authors/version only, run `julia setup.jl` |
| `Package JLD2 not found in current path` when plotting | plotting script run outside the plotting environment | `julia setup.jl plotting`, then `julia --project=plotting …` from this folder |
| `Circular dependency detected` / `Failed to precompile Makie` | Makie cannot be precompiled with Julia 1.11.0 (extension cycle) | the plotting environment uses Plots.jl instead; if an old `plotting/` folder with Makie exists, delete it and run `julia setup.jl plotting` again |
| `Assignment to prob in soft scope is ambiguous` | Julia scoping warning in scripts | harmless |
| runs much slower than expected | FFTs planned for all cores of the node | set `JULIA_CPU_THREADS` equal to `-t` (Section 4) |

---

## Appendix — batch scripts

All scripts are in `scripts/slurm/` and are submitted from the `QGTestbed` folder with `sbatch`.
Before the first use, replace `<your_account>` and `<your_qos>`, and check that `module load julia` loads the
same Julia version that was used for `julia setup.jl`. Runs started with `&` inside a script must be followed by
`wait`, otherwise the job ends immediately and Slurm stops the runs.

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
# If the job hits the time limit, submit it again: --resume skips finished runs.
#SBATCH --job-name=qg_coarse
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --time=48:00:00
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=coarse_%j.out

module load julia                  # the same Julia module used for `julia setup.jl`
module load parallel 2>/dev/null   # GNU parallel, if not on the default path
export JULIA_NUM_THREADS=1

mkdir -p logs
parallel -j 100 --joblog coarse_joblog.txt --resume --results logs/ < coarse_jobs.txt
```

### Tuning runs (Section 7b): `scripts/slurm/tune.sh`

```bash
#!/bin/bash
# Tuning runs from a job list written by scripts/tune_plan.jl, single-threaded, 100 at a time.
#   julia scripts/tune_plan.jl --closure reynolds --param C_R=1,3,7,15,30 --cases 40:64,50:64 > tune_jobs.txt
#   sbatch scripts/slurm/tune.sh                 # or: JOBS=other_list.txt sbatch scripts/slurm/tune.sh
# Resubmit to finish interrupted lists: --resume skips completed runs.
#SBATCH --job-name=qg_tune
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --time=48:00:00
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>
#SBATCH --output=tune_%j.out

module load julia                  # the same Julia module used for `julia setup.jl`
module load parallel 2>/dev/null
export JULIA_NUM_THREADS=1

JOBS=${JOBS:-tune_jobs.txt}
mkdir -p logs/tune
parallel -j 100 --joblog ${JOBS%.txt}_joblog.txt --resume --results logs/tune/ < "$JOBS"
```

---

## References

- Perezhogin, P. & Glazunov, A. (2023). Subgrid parameterizations of ocean mesoscale eddies based on Germano decomposition. JAMES. arXiv:2304.06789
- Mak, J., Marshall, D. P., Madec, G. & Maddison, J. R. (2023). Scale-awareness in an eddy energy constrained mesoscale eddy parameterization. JAMES. arXiv:2306.08988
- Jansen, M. F., Adcroft, A., Khani, S. & Kong, H. (2019). Toward an energetically consistent, resolution aware parameterization of ocean mesoscale eddies. JAMES.
- Ross, A., Li, Z., Perezhogin, P., Fernandez-Granda, C. & Zanna, L. (2023). Benchmarking of machine learning ocean subgrid parameterizations in an idealized model. JAMES.
- Constantinou, N. C. et al. (2021). GeophysicalFlows.jl: Solvers for geophysical fluid dynamics problems in periodic domains on CPUs & GPUs. JOSS 6(60), 3053.
