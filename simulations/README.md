# Simulation pipeline (false-positive control and sensitivity)

This folder contains the Python simulation pipeline behind the manuscript's false-positive-control
and detection-sensitivity results. **It is provided for transparency only.** It was run on a SLURM
compute cluster and is **not needed to regenerate any figure**: the
repository ships its cached outputs in `paper/data/`, which the `paper/code/genFig*.m` scripts read
directly.

## What it computes

Both experiments use one noise model and one SMOOTH implementation (a Python port of `SMOOTHstat`,
validated against the MATLAB code), run whole-brain: both hemispheres enter one cluster test with a
single pooled permutation null (maximum positive-cluster mass over the whole brain).

- **Noise: 1/f Gaussian random field.** Per subject, independent null noise with a power-law spectrum
  `P(λ) = (λ + median λ)^-β` (DC mode dropped) on the first 5000 graph-Laplacian eigenmodes of the
  fsaverage5 registration sphere (isotropic mesh-edge graph). Electrode covariance is
  `a·normalize(Q_e diag(P) Q_eᵀ) + (1−a)·I`, i.e. a correlated fraction `a` plus a white nugget.
  **β = 29** and **a = 0.53** come from fitting the model's electrode correlation vs. geodesic distance,
  ρ(d), to the pooled within-subject empirical ρ(d) of the primary dataset (`autocorr/fit_beta.py`,
  R² ≈ 0.97). Two nugget conditions are run: **a = 0.53** (empirical, with nugget) and **a = 1.0**
  (no nugget; pure 1/f, the adversarial, maximally autocorrelated case).
- **SMOOTH test.** Geodesic Gaussian smoothing (kernel width 12 mm), per-subject electrode-graph
  eigenbasis (σ = 8 mm), ≥ 3 subjects per vertex, cluster-forming α = 0.05, one-tailed (positive)
  α = 0.05, 1000 eigenmode sign-flip permutations.
- **Comparators.** Naive electrode shuffle (DV↔electrode permutation), and Desikan-Killiany /
  Schaefer-300 parcel-mean one-sample t-tests with a sign-flip max-statistic permutation FWER over parcels
  (same one-tailed α).
- **Specificity (FWER), `simulation/fwer.py`.** Per cell, 1000 null datasets (no effect) × 1000
  permutations; a false positive is any significant cluster/parcel. `--mode beta` sweeps 20 log-spaced
  integer β in [1, 50] (including 29) at both nugget conditions.
- **Sensitivity, `simulation/run_sensitivity.py`.** At every testable vertex (14,075 seeds) a geodesic
  Gaussian bump (FWHM 20 mm, peak SNR 4) is planted in one hemisphere on top of the a = 0.53, β = 29 noise;
  100 iterations per seed, 1000 permutations per test. Records SMOOTH sensitivity (significant cluster
  within 12 mm of the seed), off-target FPR, localization error, Dice, and the DK / Schaefer-300
  sensitivity and FPR on the same draws.
- **FDR atlas benchmark, `roi_fdr_full_sweep.py`** (with `roi_fdr_lean.py`). Same planted-effect model at
  every testable vertex, 1000 iterations per seed, atlas arms only: per-parcel one-tailed t-test with
  Benjamini-Hochberg FDR (q = 0.05) across parcels. SMOOTH values are reused from the main sensitivity run.
- **Robustness arms (Supplementary Figure 3), `robustness/`.** The published anisotropy and
  nonstationarity sweeps were run with the self-contained code in `robustness/` (see its README), at
  a = 1.0 with 1000 null datasets × 1000 permutations per value.
  *Anisotropy:* the sphere coordinates are warped by `diag(ρ^-½, ρ^½, 1)` and a separate 1500-mode
  eigenbasis is built per ρ from a Gaussian radius graph (σ = 8 mm) on the warped coordinates; fixed β = 8,
  ρ = 1, 2, 3, 5, 10.
  *Nonstationarity:* an anterior-posterior mixture of two 1/f fields with slopes β/√q (posterior, rougher)
  and β·√q (anterior, smoother) on the 1500-mode isotropic eigenbasis; β = 10, q = 1, 2, 4, 6, 8.
  Spectral exponents are not comparable across eigenbases; these values give electrode-level spatial
  autocorrelation (Moran's I ≈ 0.48 and 0.45) equivalent to the main β = 29 fields on 5000 modes (≈ 0.47).
  `fwer.py --mode aniso | nonstat` is an alternative implementation of the same two arms (per-ρ fitted β
  via `autocorr/calibrate_aniso_beta.py`; β = 29 on the 5000-mode basis for nonstationarity). It was not
  used for the figure.

## Stage → output → figure

| stage (code/ unless noted) | output | feeds |
|---|---|---|
| `autocorr/build_eigmodes.py` | `out/eigmodes/mesh_eig5000_*.npz` (large, not shipped) | all 1/f noise |
| `simulation/fwer.py --mode beta` → `simulation/merge_fwer.py` | `out/socialcog/fwer_beta/` | — |
| `viz/export_fwer_mat.py` | `paper/data/fwer_sims_1f/fwer_beta_<β>_a1.mat` (a = 1 cells) | Figure 3 |
| `simulation/run_sensitivity.py` → `simulation/merge_sensitivity.py` | `per_vertex_maps.npz` (copy in `paper/data/`) | Fig. 4, S4, S5 |
| `viz/export_benchmark_mat.py` | `paper/data/benchmark_fig_data.mat` | Figure 4 |
| `roi_fdr_full_sweep.py` → `viz/export_benchmark_fdr_mat.py` | `paper/data/benchmark_fig_data_fdr.mat` | Supp. Figure 5 |
| `viz/export_methods_fig_data.py` | `paper/data/methods_fig_data.mat` | Supp. Figure 4 |
| `robustness/run_sweeps.sh` (`robustness/fwer_sweep.py`) | `paper/data/anisononstat/fwer_sweep_*_nonug.npz` | Supp. Figure 3 |
| `viz/export_anisononstat_v3.py` | `paper/data/anisononstat_fig_data_v3.mat` | Supp. Figure 3 |
| `config.py`, `cs_common.py`, `sensitivity_core.py` | imported by `paper/code/export_fig1.py` | Figure 1 inputs |

`export_anisononstat_v3.py` reads the cached anisotropy / nonstationarity FWER sweeps in
`paper/data/anisononstat/` produced by `robustness/` (a = 1; ρ = 1, 2, 3, 5, 10 and q = 1, 2, 4, 6, 8;
1000 null datasets each; recorded at tail 'both'), reports FWER re-thresholded at the one-tailed 0.05 level, and draws the illustrative exemplar fields (its header
documents the slope and basis used for each row). `export_methods_fig_data.py` refits β and a on the active
eigenbasis (reproducing β = 29, a = 0.53) and runs one illustrative group test (300 permutations).

## Inputs

- **Electrodes / DVs:** `<repo>/demo/source.mat` (the primary dataset, 38 participants; internal label
  `socialcog`). Override with `SENS_ELEC`.
- **Surfaces:** FreeSurfer `fsaverage5` (`lh/rh.pial`, `.sphere.reg`, `.inflated`, `.curv`), read from
  `$SUBJECTS_DIR/fsaverage5/surf` or `$FREESURFER_HOME/subjects/fsaverage5/surf` (override: `SENS_FS_SURF`).
- **Atlas lookups:** `<repo>/paper/data/{aparc,schaefer300}_fs5.npz` (per-vertex parcel labels).

All parameters live in `code/config.py` (locked to the values above, each env-overridable). Outputs are
written under `simulations/out/`; exporters write to `paper/data/` unless their `SENS_*_OUT` override is set.

## How it was run

Python 3.11 with `requirements.txt` (numpy, scipy, matplotlib), inside the container defined in
`env/apptainer.def` (or `env/Dockerfile`); the two exporters that read curvature also need `nibabel`.

```bash
cd simulations
apptainer build env/smooth-simulations.sif env/apptainer.def   # optional; else set PYTHON_SETUP
export FREESURFER_HOME=<freesurfer>  CPUS=<cores/task>  PARTITION=<partition>  ACCOUNT=<account>
./slurm/submit_all.sh            # (0) eigenmodes -> (1) FWER sweeps + merge -> (2) sensitivity + merge
SENS_MESH_EIG_TAG=5000 python code/autocorr/calibrate_aniso_beta.py   # only for the alternative fwer.py aniso mode
NCORES=<cores> ./robustness/run_sweeps.sh                             # Supp. Fig. 3 sweeps (one multi-core machine)
python code/roi_fdr_full_sweep.py <nworkers>                          # FDR sweep (one multi-core machine)
python code/viz/export_*.py                                           # cached .mat files for paper/data/
```

`slurm/env.sh` holds every site-specific setting (Python activation, partition, account, cores, memory,
wall time, container image); array sizes are computed from the parameter grids. `submit_all.sh eig | fwer
| sens` runs a single stage. Each array task is resumable and skips completed outputs.
