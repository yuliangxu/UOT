#!/usr/bin/env bash
#SBATCH --job-name=uot-clt
#SBATCH --array=1-100
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=/cwork/yx306/UOT/out/uot-clt-%A_%a.out
#SBATCH --error=/cwork/yx306/UOT/out/uot-clt-%A_%a.err
#SBATCH --partition=mastatlab
#SBATCH --account=mastatlab

set -euo pipefail

: "${SLURM_ARRAY_TASK_ID:?This file must run as a SLURM array job}"

REPO_ROOT="${REPO_ROOT:-/hpc/home/yx306/UOT}"
VENV_PYTHON="${VENV_PYTHON:-/hpc/home/yx306/venvs/uot/bin/python}"
CLT_RUN_ROOT="${CLT_RUN_ROOT:-/cwork/yx306/UOT/clt-gaussian-10000}"
CLT_SHARD_DIR="${CLT_SHARD_DIR:-${CLT_RUN_ROOT}/shards}"
CLT_REPS="${CLT_REPS:-10000}"
CLT_REP_START="${CLT_REP_START:-1}"
CLT_NUM_SHARDS="${CLT_NUM_SHARDS:-100}"

if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > CLT_NUM_SHARDS )); then
  echo "SLURM_ARRAY_TASK_ID must be between 1 and CLT_NUM_SHARDS" >&2
  exit 2
fi

cd "$REPO_ROOT"
mkdir -p "$CLT_SHARD_DIR"

export RETICULATE_PYTHON="$VENV_PYTHON"
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "Repository: $REPO_ROOT"
echo "Python: $RETICULATE_PYTHON"
echo "Replications: $CLT_REP_START:$CLT_REPS"
echo "Shard: $SLURM_ARRAY_TASK_ID/$CLT_NUM_SHARDS"

Rscript --vanilla simulation/ET_round2/CLT_gaussian_hpc.R \
  --mode=shard \
  --reps="$CLT_REPS" \
  --rep-start="$CLT_REP_START" \
  --num-shards="$CLT_NUM_SHARDS" \
  --shard-id="$SLURM_ARRAY_TASK_ID" \
  --shard-dir="$CLT_SHARD_DIR"
