#!/usr/bin/env bash
#SBATCH --job-name=sim-sens
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=logs/sens_%A_%a.out
# Whole-brain sensitivity, one array task = one seed chunk (strided). Output label via SENS_DATASET
# (default socialcog), number of chunks via NCHUNKS (array size).
#   SENS_DATASET=socialcog NCHUNKS=200 sbatch --array=0-199 --cpus-per-task=$CPUS --mem=$MEM --time=$TIME slurm/sensitivity.sbatch
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/env.sh"
mkdir -p "$PROJ/logs"
RUN "$PROJ/code/simulation/run_sensitivity.py" \
    --chunk "${SLURM_ARRAY_TASK_ID:-0}" --nchunks "${NCHUNKS:-200}" \
    --ncores "${SLURM_CPUS_PER_TASK:-$CPUS}" --niter "${SENS_NITER:-100}" --nperm "${SENS_NPERM:-1000}"
