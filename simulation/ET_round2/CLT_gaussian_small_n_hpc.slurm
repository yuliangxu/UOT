#!/usr/bin/env bash
#SBATCH --job-name=uot-clt-small-n
#SBATCH --array=1-100
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/cwork/yx306/UOT/out/uot-clt-small-n-%A_%a.out
#SBATCH --error=/cwork/yx306/UOT/out/uot-clt-small-n-%A_%a.err
#SBATCH --partition=mastatlab
#SBATCH --account=mastatlab

set -euo pipefail

: "${SLURM_ARRAY_TASK_ID:?This file must run as a SLURM array job}"

REPO_ROOT="${REPO_ROOT:-/hpc/home/yx306/OT/UOT}"
R_BIN="${R_BIN:-/opt/apps/rhel9/R-4.4.3/bin}"
VENV_PYTHON="${VENV_PYTHON:-/hpc/home/yx306/venvs/uot/bin/python}"
CLT_RUN_ROOT="${CLT_RUN_ROOT:-/cwork/yx306/UOT/out/clt-gaussian-10000}"
CLT_SMALL_N_SHARD_DIR="${CLT_SMALL_N_SHARD_DIR:-${CLT_RUN_ROOT}/small-n-shards}"
CLT_REPS="${CLT_REPS:-10000}"
CLT_REP_START="${CLT_REP_START:-1}"
CLT_NUM_SHARDS="${CLT_NUM_SHARDS:-100}"

if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > CLT_NUM_SHARDS )); then
  echo "SLURM_ARRAY_TASK_ID must be between 1 and CLT_NUM_SHARDS" >&2
  exit 2
fi

case "$CLT_SMALL_N_SHARD_DIR" in
  /cwork/yx306/UOT/out/*) ;;
  *)
    echo "CLT_SMALL_N_SHARD_DIR must be under /cwork/yx306/UOT/out" >&2
    exit 2
    ;;
esac

cd "$REPO_ROOT"
mkdir -p "$CLT_SMALL_N_SHARD_DIR"

export PATH="$R_BIN:$PATH"
export RETICULATE_PYTHON="$VENV_PYTHON"
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "Repository: $REPO_ROOT"
echo "Sample sizes: 10,20"
echo "Replications: $CLT_REP_START:$CLT_REPS"
echo "Shard: $SLURM_ARRAY_TASK_ID/$CLT_NUM_SHARDS"

Rscript --vanilla simulation/ET_round2/CLT_gaussian_small_n_hpc.R \
  --mode=shard \
  --reps="$CLT_REPS" \
  --rep-start="$CLT_REP_START" \
  --num-shards="$CLT_NUM_SHARDS" \
  --shard-id="$SLURM_ARRAY_TASK_ID" \
  --shard-dir="$CLT_SMALL_N_SHARD_DIR"
