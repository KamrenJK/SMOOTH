#!/usr/bin/env bash
# One-shot submit of the whole pipeline with correct dependencies. Resources come from env.sh (all
# overridable). Stages: (0) build shared eigenmodes -> (1) FWER sweeps beta/aniso/nonstat -> (2)
# whole-brain sensitivity -> (3) merge+plot. Nothing is hard-coded to a node count; array
# sizes are computed from the actual grids, so this scales to whatever allocation you have.
#
#   ./slurm/submit_all.sh                 # everything
#   ./slurm/submit_all.sh fwer            # just the FWER sweeps (assumes eigenmodes already built)
#   ./slurm/submit_all.sh sens            # just sensitivity
# Env knobs (see slurm/env.sh): CPUS MEM TIME PARTITION ACCOUNT IMAGE PYTHON_SETUP SENS_FS_SURF ;
#                                SENS_NSIM SENS_NPERM SENS_NITER
#                                SENS_FWER_CHUNKS (split each FWER cell across more tasks) ; NCHUNKS (sens)
set -euo pipefail
cd "$(dirname "$0")/.."
source slurm/env.sh
STAGE="${1:-all}"
FLAGS=$(sbatch_flags)
COMMON=(--cpus-per-task "$CPUS" --mem "$MEM" --time "$TIME" $FLAGS)
mkdir -p logs

submit() { sbatch --parsable "$@" | tail -n1; }   # echo the job id

EIG_DEP=""
if [[ "$STAGE" == "all" || "$STAGE" == "eig" ]]; then
  echo "[0] build_eigmodes (iso + aniso rho=1..10)"
  EIG=$(submit --array=0-10 "${COMMON[@]}" slurm/build_eigmodes.sbatch)
  echo "    job $EIG"; EIG_DEP="--dependency=afterok:$EIG"
fi

if [[ "$STAGE" == "all" || "$STAGE" == "fwer" ]]; then
  for MODE in beta aniso nonstat; do
    N=$(SENS_FWER_MODE=$MODE RUN code/simulation/fwer.py --mode "$MODE" --nugs "${SENS_NUGS:-0.53,1.0}" --count)
    echo "[1] FWER $MODE (nugs=${SENS_NUGS:-0.53,1.0}) : $N array tasks"
    JID=$(SENS_FWER_MODE=$MODE SENS_NUGS="${SENS_NUGS:-0.53,1.0}" submit --array=0-$((N-1)) $EIG_DEP "${COMMON[@]}" slurm/fwer.sbatch)
    echo "    job $JID"
    sbatch --parsable --dependency=afterok:$JID $FLAGS --time=00:20:00 --cpus-per-task=2 --mem=8G \
      --wrap "source slurm/env.sh; RUN code/simulation/merge_fwer.py --mode $MODE" >/dev/null
  done
fi

if [[ "$STAGE" == "all" || "$STAGE" == "sens" ]]; then
  NCH="${NCHUNKS:-200}"
  for DS in "${SENS_DATASET:-socialcog}"; do
    echo "[2] sensitivity $DS : $NCH chunks"
    JID=$(SENS_DATASET=$DS NCHUNKS=$NCH submit --array=0-$((NCH-1)) $EIG_DEP "${COMMON[@]}" slurm/sensitivity.sbatch)
    echo "    job $JID"
    SENS_DATASET=$DS sbatch --parsable --dependency=afterok:$JID $FLAGS --time=00:30:00 --cpus-per-task=2 --mem=16G \
      --wrap "source slurm/env.sh; SENS_DATASET=$DS RUN code/simulation/merge_sensitivity.py" >/dev/null
  done
fi
echo "submitted ($STAGE)."
