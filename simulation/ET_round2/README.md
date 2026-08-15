# ET round 2: Gaussian UOT CLT experiment

This folder contains the numerical experiment for Theorem 1 and a
scheduler-safe runner for a larger Monte Carlo study.

- `CLT_gaussian_experiment.R` runs or resumes the experiment serially and
  creates the result tables and plots.
- `CLT_gaussian_hpc.R` runs independent immutable shards, validates and
  merges them, and is safe to use with a scheduler array.
- `CLT_gaussian_experiment.md` describes the experiment and interprets the
  current plots.
- `CLT_gaussian_hpc.slurm` is a SLURM array-job template.

Generated `result/` and `plot/` directories are intentionally ignored by
Git. Run every command below from the repository root.

## Software environment

The R requirements are `MASS`, `ggplot2`, and `reticulate`. The selected
Python must contain POT (imported as `ot`), NumPy, and SciPy. On a cluster,
create the Python environment on a login or build node rather than inside
each array task. For example:

```bash
python3 -m venv /path/to/venvs/uot
/path/to/venvs/uot/bin/python -m pip install --upgrade pip
/path/to/venvs/uot/bin/python -m pip install \
  'POT==0.9.7.post1' 'numpy==2.0.2' 'scipy==1.13.1'
```

Install the R packages once in the R library used by the jobs:

```r
install.packages(c("MASS", "ggplot2", "reticulate"))
```

Set the absolute Python path before starting R. Setting numerical-library
thread counts to one prevents every array task from creating additional
BLAS/OpenMP threads.

```bash
export RETICULATE_PYTHON=/path/to/venvs/uot/bin/python
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

Rscript --vanilla simulation/ET_round2/CLT_gaussian_hpc.R --mode=preflight
```

The preflight imports the required software, runs a small POT solve, and
checks the Gaussian population transported mass against the reference value
`0.588149213661207` for $(\epsilon,\rho)=(0.5,1)$.

## Production design

The production target is 10,000 evaluation replications at each
$n\in\{50,100,200,500\}$ and a separate 10,000-replication reference batch
at $n=500$. Thus `--reps=10000` means 50,000 UOT fits in total, not 10,000.
The reference and evaluation batches use disjoint deterministic seed streams.

The shard runner preserves the global replication ID in every job. Results
therefore do not depend on the number of shards, retry order, or scheduler.
Each task writes only its own complete shard file; no two tasks write a
shared checkpoint.

### Fresh run

First test the workflow in scratch with a tiny design:

```bash
Rscript --vanilla simulation/ET_round2/CLT_gaussian_hpc.R \
  --mode=shard --smoke-test --reps=4 --rep-start=1 \
  --num-shards=2 --shard-id=1 \
  --shard-dir=/scratch/$USER/uot_clt_smoke/shards

Rscript --vanilla simulation/ET_round2/CLT_gaussian_hpc.R \
  --mode=shard --smoke-test --reps=4 --rep-start=1 \
  --num-shards=2 --shard-id=2 \
  --shard-dir=/scratch/$USER/uot_clt_smoke/shards

Rscript --vanilla simulation/ET_round2/CLT_gaussian_hpc.R \
  --mode=merge --smoke-test --reps=4 --rep-start=1 \
  --num-shards=2 \
  --shard-dir=/scratch/$USER/uot_clt_smoke/shards \
  --output-dir=/scratch/$USER/uot_clt_smoke/result \
  --plot-dir=/scratch/$USER/uot_clt_smoke/plot
```

For a fresh production run with 100 shards, submit array indices 1--100.
The supplied SLURM template reads its paths and design from environment
variables:

```bash
mkdir -p /cwork/yx306/UOT/out

# These defaults are already in the template; override them if needed.
export REPO_ROOT=/hpc/home/yx306/UOT
export VENV_PYTHON=/hpc/home/yx306/venvs/uot/bin/python
export CLT_RUN_ROOT=/cwork/yx306/UOT/clt-gaussian-10000
export CLT_SHARD_DIR="$CLT_RUN_ROOT/shards"
export CLT_RESULT_DIR="$CLT_RUN_ROOT/result"
export CLT_PLOT_DIR="$CLT_RUN_ROOT/plot"
export CLT_REPS=10000
export CLT_REP_START=1
export CLT_NUM_SHARDS=100

Rscript --vanilla simulation/ET_round2/CLT_gaussian_hpc.R \
  --mode=preflight --shard-dir="$CLT_SHARD_DIR" \
  --output-dir="$CLT_RESULT_DIR" --plot-dir="$CLT_PLOT_DIR"

sbatch --array=1-100 simulation/ET_round2/CLT_gaussian_hpc.slurm
```

The template uses the `mastatlab` partition/account, 8 GB of memory, one CPU,
and no GPU. POT is CPU-based here, and each task explicitly limits BLAS and
OpenMP to one thread. Adjust the 24-hour limit after timing the smoke run if
your cluster requires a different request. The `/cwork/yx306/UOT/out`
directory must exist before `sbatch`, because SLURM opens its log files before
the script starts.

After every array task succeeds, merge once. Merge validates all inputs and
then automatically invokes the ordinary analysis to create metadata, the
covariance object, CSV summaries, and PDFs:

```bash
Rscript --vanilla simulation/ET_round2/CLT_gaussian_hpc.R \
  --mode=merge --reps=10000 --rep-start=1 --num-shards=100 \
  --shard-dir="$CLT_SHARD_DIR" --output-dir="$CLT_RESULT_DIR" \
  --plot-dir="$CLT_PLOT_DIR"
```

Merge mode fails without modifying the canonical checkpoints if a shard is
missing, incompatible, nonconverged, or has an incorrect job/seed set. Keep
the shard files and the generated merge manifest until the analysis is
archived.

### Reuse the completed 2,000 replications

To compute only replications 2,001--10,000, first copy these four compatible
files from the existing local `result/` directory into `CLT_RESULT_DIR`:

```text
clt_gaussian_all_projections_checkpoint.rds
clt_gaussian_all_projections_independent_reference_checkpoint.rds
clt_gaussian_all_projections_metadata.rds
clt_gaussian_all_projections_independent_reference_metadata.rds
```

Back them up and record checksums before transfer. Then use
`CLT_REP_START=2001` in the array job and use `--rep-start=2001` in merge
mode. Keep the default master seed `20260814`. The merge requires a complete
base for replication IDs 1--2,000 and rejects missing or mismatched jobs.

Do not run the serial script or multiple array tasks against the same shard
ID at the same time. A retried task is safe: it validates and reuses its
already-complete immutable shard.

## Completion checks

For the production design, the merged evaluation checkpoint must contain
40,000 rows and the independent reference checkpoint 10,000 rows:

```r
evaluation <- readRDS(file.path(
  Sys.getenv("CLT_RESULT_DIR"),
  "clt_gaussian_all_projections_checkpoint.rds"
))
reference <- readRDS(file.path(
  Sys.getenv("CLT_RESULT_DIR"),
  "clt_gaussian_all_projections_independent_reference_checkpoint.rds"
))

stopifnot(
  nrow(evaluation$diagnostics) == 40000L,
  all(table(evaluation$diagnostics$n) == 10000L),
  nrow(reference$diagnostics) == 10000L,
  identical(unique(reference$diagnostics$n), 500L),
  all(evaluation$diagnostics$converged),
  all(reference$diagnostics$converged),
  !anyDuplicated(evaluation$diagnostics$seed),
  !anyDuplicated(reference$diagnostics$seed),
  length(intersect(
    evaluation$diagnostics$seed,
    reference$diagnostics$seed
  )) == 0L
)
```

Archive both checkpoints and metadata files, the merge manifest, covariance
RDS, CSV summaries, all generated PDFs, scheduler logs, and the recorded
R/Python/package versions.

Increasing from 2,000 to 10,000 replications reduces Monte Carlo covariance
estimation error by roughly $\sqrt{2000/10000}=1/\sqrt{5}$. It does not
increase $n$, make the $n=500$ reference covariance analytic, or prove
convergence to the unknown limiting covariance. Interpret it as a more
precise finite-$n$ Monte Carlo comparison.
