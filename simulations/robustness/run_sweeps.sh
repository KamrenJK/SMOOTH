#!/usr/bin/env bash
# Exact commands for the two FWER sweeps behind Supplementary Figure 3 (anisotropic and nonstationary
# 1/f noise, no nugget). Outputs: out/fwer_sweep_aniso1f_nonug.npz and out/fwer_sweep_nonstat1f_nonug.npz
# (cached copies in <repo>/paper/data/anisononstat/). Each sweep is 5 values x 1000 null datasets x 1000
# permutations; it was run on a 32-core node.
#   ./run_sweeps.sh            # cache eigenmodes, then both sweeps
#   NCORES=64 ./run_sweeps.sh  # worker count (does not affect results; trials are seeded per (value, trial))
# Surfaces: fsaverage5 from $SENS_FS_SURF, $SUBJECTS_DIR or $FREESURFER_HOME (see config.py).
# Electrodes: <repo>/demo/source.mat (override SENS_ELEC).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"                               # simulations/robustness/

# Optional Python environment hook (numpy + scipy), e.g. PYTHON_SETUP='source <venv>/bin/activate'
if [[ -n "${PYTHON_SETUP:-}" ]]; then set +u; eval "$PYTHON_SETUP"; set -u; fi

NCORES="${NCORES:-32}"
export SENS_DATASET=socialcog                                     # output label recorded in the sweep files
export SENS_CLUSTERALPHA=0.05                                     # cluster-forming alpha
export SENS_TAIL=both                                             # FWER is re-thresholded one-tailed at 0.05 by the exporter
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1   # 1 BLAS thread/worker

# (1) eigenmode caches (idempotent; all cores for the one-time eigsh)
OMP_NUM_THREADS=$NCORES OPENBLAS_NUM_THREADS=$NCORES MKL_NUM_THREADS=$NCORES python cache_mesh_eigmodes.py 1500
OMP_NUM_THREADS=$NCORES OPENBLAS_NUM_THREADS=$NCORES MKL_NUM_THREADS=$NCORES python cache_mesh_eigmodes_aniso.py 1500 8 1 2 3 5 10

# (2) the two published sweeps
python fwer_sweep.py --noise aniso1f --rhos 1 2 3 5 10 --a 1.0 --beta 8 --outtag aniso1f_nonug --nsim 1000 --nperm 1000 --ncores "$NCORES"
python fwer_sweep.py --noise nonstat1f --qs 1 2 4 6 8 --a 1.0 --beta 10 --outtag nonstat1f_nonug --nsim 1000 --nperm 1000 --ncores "$NCORES"
