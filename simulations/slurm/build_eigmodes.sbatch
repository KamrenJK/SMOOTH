#!/usr/bin/env bash
#SBATCH --job-name=sim-eig
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=logs/eig_%A_%a.out
# Setup stage: build the shared fs5 graph-Laplacian eigenmode caches (isotropic + warped rho=1..10)
# into $PROJ/out/eigmodes. Dataset-independent -> run ONCE. Array index 0=iso, 1..10=aniso rho.
# Submit with:  sbatch --array=0-10 --cpus-per-task=$CPUS --mem=$MEM --time=$TIME slurm/build_eigmodes.sbatch
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/env.sh"
mkdir -p "$PROJ/logs"
# The eigendecomposition is a SINGLE process (scipy eigh) -- override env.sh's 1-thread pin so BLAS
# uses the whole node. (The sims keep 1 thread/worker; this override is local to the eig build.)
NT="${SLURM_CPUS_PER_TASK:-$CPUS}"
export OMP_NUM_THREADS="$NT" OPENBLAS_NUM_THREADS="$NT" MKL_NUM_THREADS="$NT" NUMEXPR_NUM_THREADS="$NT"
i="${SLURM_ARRAY_TASK_ID:-0}"
if [[ "$i" == "0" ]]; then
  RUN "$PROJ/code/autocorr/build_eigmodes.py" iso
else
  RUN "$PROJ/code/autocorr/build_eigmodes.py" aniso "$i"
fi
