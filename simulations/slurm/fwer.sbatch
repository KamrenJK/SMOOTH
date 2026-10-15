#!/usr/bin/env bash
#SBATCH --job-name=sim-fwer
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=logs/fwer_%A_%a.out
# FWER sweep, one array task = one (noise-cell x nugget x sim-chunk). MODE (beta|aniso|nonstat) is passed
# via the SENS_FWER_MODE env. Array size = Ncells*CHUNKS (submit_all.sh computes it).
#   SENS_FWER_MODE=beta sbatch --array=0-39 --cpus-per-task=$CPUS --mem=$MEM --time=$TIME slurm/fwer.sbatch
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/env.sh"
mkdir -p "$PROJ/logs"
MODE="${SENS_FWER_MODE:-beta}"
RUN "$PROJ/code/simulation/fwer.py" --mode "$MODE" --task "${SLURM_ARRAY_TASK_ID:-all}" \
    --nugs "${SENS_NUGS:-0.53,1.0}" \
    --nsim "${SENS_NSIM:-1000}" --nperm "${SENS_NPERM:-1000}" --ncores "${SLURM_CPUS_PER_TASK:-$CPUS}"
