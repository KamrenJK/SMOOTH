#!/usr/bin/env bash
# Shared cluster config for the simulation pipeline, sourced by every sbatch script and submit_all.sh.
# Everything is env-overridable so the same scripts run on any allocation.
set -euo pipefail

# --- paths: config.py derives everything from these ---
export PROJ="${PROJ:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"        # simulations/
export SMOOTH_REPO="${SMOOTH_REPO:-$(cd "$PROJ/.." && pwd)}"                     # repository root
export PYTHONPATH="$PROJ/code:${PYTHONPATH:-}"
# fsaverage5 surf/ directory of a FreeSurfer installation (lh/rh .pial, .sphere.reg, .inflated, .curv)
if [[ -z "${SENS_FS_SURF:-}" ]]; then
  if [[ -n "${SUBJECTS_DIR:-}" && -d "$SUBJECTS_DIR/fsaverage5/surf" ]]; then
    export SENS_FS_SURF="$SUBJECTS_DIR/fsaverage5/surf"
  else
    export SENS_FS_SURF="${FREESURFER_HOME:-/usr/local/freesurfer}/subjects/fsaverage5/surf"
  fi
fi

# --- container: set IMAGE to a .sif to run inside Apptainer; leave empty to use the host python ---
export IMAGE="${IMAGE:-$PROJ/env/smooth-simulations.sif}"

# --- Python environment (numpy/scipy/matplotlib) ------------------------------------------------
# Batch jobs start in a FRESH shell, so whatever makes `python` work interactively (conda/module/venv)
# must be re-activated here -- this file is sourced on BOTH the login node (submit_all --count) and the
# compute nodes (the jobs). Skipped when the Apptainer image ($IMAGE) exists (the container ships them).
# Set PYTHON_SETUP to the activation command for your cluster, e.g.
#     PYTHON_SETUP='module load python'
#     PYTHON_SETUP='source <conda-root>/etc/profile.d/conda.sh && conda activate <env>'
#     PYTHON_SETUP='source <venv>/bin/activate'
if [[ ! -f "$IMAGE" && -n "${PYTHON_SETUP:-}" ]]; then
  set +u                                                     # module/conda init aren't `set -u` clean
  eval "$PYTHON_SETUP"
  set -u
fi

# --- resources (flexible; SLURM sets SLURM_CPUS_PER_TASK, which config.py reads as NCORES) ---
export PARTITION="${PARTITION:-}"          # <partition> (blank -> cluster default)
export CPUS="${CPUS:-128}"                 # cores per array task
export MEM="${MEM:-64G}"                   # RAM per task (5000-mode caches ~0.4 GB/worker)
export TIME="${TIME:-08:00:00}"
export ACCOUNT="${ACCOUNT:-}"              # <account> (blank -> none)

# --- pin BLAS to 1 thread/worker so N workers use N cores (config also sets these) ---
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg

# RUN <python-args...> : exec the driver, inside the container if IMAGE exists, else host python.
# The container sees the repository and the fsaverage5 surfaces at their host paths.
RUN() {
  if [[ -f "$IMAGE" ]]; then
    apptainer exec --bind "$SMOOTH_REPO" --bind "$SENS_FS_SURF" \
      --env "PROJ=$PROJ,SMOOTH_REPO=$SMOOTH_REPO,SENS_FS_SURF=$SENS_FS_SURF,PYTHONPATH=$PROJ/code" \
      "$IMAGE" python "$@"
  else
    python "$@"
  fi
}
export -f RUN

# SBATCH_FLAGS: assemble optional sbatch flags (partition/account) only when set.
sbatch_flags() {
  local f=()
  [[ -n "$PARTITION" ]] && f+=(--partition "$PARTITION")
  [[ -n "$ACCOUNT"   ]] && f+=(--account "$ACCOUNT")
  echo "${f[@]:-}"
}
