# Robustness sweeps (Supplementary Figure 3)

This folder holds the code that produced the two cached FWER sweeps behind Supplementary Figure 3
(false-positive control under anisotropic and nonstationary noise):

- `paper/data/anisononstat/fwer_sweep_aniso1f_nonug.npz`
- `paper/data/anisononstat/fwer_sweep_nonstat1f_nonug.npz`

`simulations/code/viz/export_anisononstat_v3.py` reads them to build `paper/data/anisononstat_fig_data_v3.mat`,
which `paper/code/genFigS3.m` draws. This folder is self-contained: `fwer_sweep.py` (driver),
`sensitivity_core.py` + `config.py` (engine), `cache_mesh_eigmodes.py` / `cache_mesh_eigmodes_aniso.py`
(eigenbases), and `run_sweeps.sh` (the exact commands). It is separate from `simulations/code/`.

## Settings

- No nugget (a = 1.0). 1,000 null datasets × 1,000 permutations per swept value.
- SMOOTH test: kernel width 12 mm, graph sigma 8 mm, cluster-forming alpha 0.05, both hemispheres in one
  test with a pooled null. Significance masks were recorded at tail `both`; the exporter re-thresholds FWER
  one-tailed at 0.05.
- **Anisotropy:** a separate 1,500-mode eigenbasis per ρ, from a Gaussian-weighted radius graph
  (sigma 8 mm) on sphere coordinates warped by `diag(ρ^-½, ρ^½, 1)`; β = 8; ρ = 1, 2, 3, 5, 10.
- **Nonstationarity:** an anterior-posterior mixture of two 1/f fields with slopes β/√q and β·√q on the
  1,500-mode isotropic eigenbasis; β = 10; q = 1, 2, 4, 6, 8.

Spectral exponents are not comparable across eigenbases; these values give spatial autocorrelation
equivalent to the main simulations (electrode-level Moran's I ≈ 0.48 and 0.45, vs 0.47 for the main
β = 29 fields on the 5,000-mode basis).

## Running

```bash
cd simulations/robustness
export FREESURFER_HOME=<freesurfer>      # fsaverage5 surfaces (or SUBJECTS_DIR / SENS_FS_SURF)
NCORES=<cores> ./run_sweeps.sh           # optional PYTHON_SETUP='<activation command>'
```

`run_sweeps.sh` sets `SENS_CLUSTERALPHA=0.05`, `SENS_TAIL=both` and single-threaded BLAS, builds the
eigenmode caches (`cache_mesh_eigmodes.py 1500`; `cache_mesh_eigmodes_aniso.py 1500 8 1 2 3 5 10`), then runs:

```bash
python fwer_sweep.py --noise aniso1f --rhos 1 2 3 5 10 --a 1.0 --beta 8 --outtag aniso1f_nonug --nsim 1000 --nperm 1000
python fwer_sweep.py --noise nonstat1f --qs 1 2 4 6 8 --a 1.0 --beta 10 --outtag nonstat1f_nonug --nsim 1000 --nperm 1000
```

Electrodes are read from `<repo>/demo/source.mat` (override `SENS_ELEC`). Eigenmode caches and sweep
results are written to `out/` (override `SENS_OUT`). Needs Python 3 with numpy and scipy.
