# ET round 2: Gaussian UOT CLT experiment

This folder contains the numerical experiment for Theorem 1 and a
scheduler-safe runner for a larger Monte Carlo study.

- `CLT_gaussian_experiment.R` runs or resumes the experiment serially and
  creates the result tables and plots.
- `CLT_gaussian_hpc.R` runs independent immutable shards, validates and
  merges them, and is safe to use with a scheduler array.
- `CLT_gaussian_summary.R` validates completed shards and creates compact
  CSV/RDS summaries, a Markdown report, and diagnostic figures without
  rerunning the simulation.
- `CLT_gaussian_experiment.md` describes the experiment and interprets the
  current plots.
- `CLT_gaussian_hpc.slurm` is a SLURM array-job template.
- `CLT_gaussian_small_n_hpc.R` and `.slurm` add evaluation-only runs at
  `n=10,20` without recomputing the existing larger sample sizes or reference.
- `CLT_gaussian_small_n_summary.R` combines both immutable runs and creates a
  relative Frobenius-distance plot, five distributional diagnostics, and a
  Jarque–Bera table. The independent `n=500` draws supply the reference
  covariance and normal overlay scales.

Generated `result/` and `plot/` directories are intentionally ignored by
Git. Run every command below from the repository root.

## Reproduce the six-size Gaussian diagnostics

Use the repository root as the working directory. Set `CLT_RUN_ROOT` to the
directory containing the completed `shards/` and `small-n-shards/` trees.
Set `CLT_ARTIFACT_DIR` to the destination for the new plots and tables.
With R and ggplot2 installed, this one command reads the base and
continuation runs and generates the six files linked below.

```bash
export CLT_RUN_ROOT=/path/to/completed/gaussian-run
export CLT_ARTIFACT_DIR="$PWD/result/gaussian-replay"

Rscript --vanilla simulation/ET_round2/CLT_gaussian_small_n_summary.R \
  --reps=100000 --base-reps=10000 --extension-rep-start=10001 \
  --base-num-shards=100 --extension-num-shards=100 \
  --base-shard-dir="$CLT_RUN_ROOT/shards" \
  --small-n-shard-dir="$CLT_RUN_ROOT/small-n-shards" \
  --output-dir="$CLT_ARTIFACT_DIR" --plot-dir="$CLT_ARTIFACT_DIR" \
  --mc-block-size=2000
```

Each evaluation uses $n$ observations **per group**, for
$n=10,20,50,100,200,500$. All five PDFs show the six sample sizes.
For every test function $f$, the fluctuation is
$\sqrt n\,\{\int f\,d\widehat\gamma_n-\int f\,d\gamma_0\}$.
The fits use the raw transport plan, squared Euclidean cost, $\epsilon=0.5$,
$\rho=1$, and KL regularization relative to the product of the marginals.
The default master seed is `20260814`. Evaluation uses stream 0 and the
independent $n=500$ reference uses stream 1.

| Output | Computation |
| --- | --- |
| [Representative fluctuations](gaussian_diagnostics/clt_gaussian_all_projections_R100000_representative_fluctuations.pdf) | Evaluation replication ID 1 at each $n$. Each heatmap shows the $7\times7$ bump centers for one of p11, p12, p21, and p22. |
| [Mean fluctuations](gaussian_diagnostics/clt_gaussian_all_projections_R100000_mean_fluctuations.pdf) | Mean of each of the 196 bump fluctuations across 100,000 evaluation replications at each $n$. |
| [SD of fluctuations](gaussian_diagnostics/clt_gaussian_all_projections_R100000_sd_fluctuations.pdf) | Sample standard deviation of each bump fluctuation across the same evaluation replications. |
| [Histograms](gaussian_diagnostics/clt_gaussian_all_projections_R100000_histograms.pdf) | 25 bins per panel for total transported mass and the central bump of each projection. Each zero-mean normal curve uses that feature's sample SD in the 100,000 independent reference draws. |
| [Normal Q–Q plots](gaussian_diagnostics/clt_gaussian_all_projections_R100000_qqplots.pdf) | The same five evaluation fluctuations, divided by their respective reference SDs and sorted within each panel, versus `qnorm(ppoints(100000))`. |
| [Jarque–Bera table](gaussian_diagnostics/jarque_bera_first3000.csv) | One row per sample size and feature, using evaluation IDs 1 through 3,000 at each size (18,000 fits). |

For the table, $m_k=R^{-1}\sum_i(x_i-\bar x)^k$ with $R=3000$.
The statistic is $R\{(m_3/m_2^{3/2})^2+(m_4/m_2^2-3)^2/4\}/6$.
The raw p-value is the upper tail of $\chi^2_2$. A single Holm correction
covers all 30 tests in increasing $n$ and the feature order `total_mass`,
`p11_bump_s04_t04`, `p12_bump_s04_t04`, `p21_bump_s04_t04`,
`p22_bump_s04_t04`. The CSV retains 17 significant digits. For the
production seed, the summary checks the JB statistic and Holm p-value in
every cell against the published rounding.

The summary validates 100,000 converged fits at each evaluation size
(600,000 total) and 100,000 independent reference fits. It checks
replication IDs, uniqueness and separation of the seed streams, finite
integrals, and solver convergence. The previous six-size summary, report,
covariance-distance PDF/CSV, and calibration RDS are also generated.

### Production verification

The archived run retained its PDFs and numerical summaries. We regenerated its
missing shards with the original workers and seed streams, then rendered all
five published PDFs and the CSV from a fresh checkout. All six evaluation sizes
had 100,000 successful fits each. The independent reference had 100,000 at
$n=500$. All 700,000 seeds were distinct and every fit had 197 finite integrals.
The six covariance distances and Monte Carlo classifications matched the archive.

A separate calculation from evaluation IDs 1–3,000 reproduced all 30 published
JB statistics and all 30 Holm-adjusted p-values at their displayed precision.
The Q–Q PDF contains all 100,000 order statistics in each of its 30 panels.
A second end-to-end run on the same six sample sizes produced the five
`R4` PDFs, `jarque_bera_first4.csv`, and the original covariance-distance plot.

The verified fit environment used R 4.4.3, MASS 7.3-65, ggplot2 4.0.2,
reticulate 1.45.0, Python 3.9, and the POT, NumPy, and SciPy versions below.

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
Rscript --vanilla simulation/ET_round2/CLT_gaussian_hpc.R --mode=preflight \
  --shard-dir="$CLT_RUN_ROOT/shards" --output-dir="$CLT_RUN_ROOT/result" \
  --plot-dir="$CLT_RUN_ROOT/plot"
```

The preflight imports the required software, runs a small POT solve, and
checks the Gaussian population transported mass against the reference value
`0.588149213661207` for $(\epsilon,\rho)=(0.5,1)$.

## Generate missing simulation inputs

Use the software environment above and a run directory accessible to all
workers. Each worker can validate and reuse a completed shard. The first
worker fits $n=50,100,200,500$ and the independent reference. The small-$n$
worker fits $n=10,20$. Each worker derives its seed from the master seed,
sample size, replication ID and stream.
Set `PARTITION` and `ACCOUNT` to your Slurm partition and account. Submit
from the repository root.

```bash
export REPO_ROOT="$PWD"
export R_BIN="$(dirname "$(command -v Rscript)")"
export VENV_PYTHON="$RETICULATE_PYTHON"
export CLT_SHARD_DIR="$CLT_RUN_ROOT/shards"
export CLT_SMALL_N_SHARD_DIR="$CLT_RUN_ROOT/small-n-shards"
export CLT_NUM_SHARDS=100
mkdir -p "$CLT_RUN_ROOT/logs"

for CLT_REP_START in 1 10001; do
  export CLT_REP_START
  if [ "$CLT_REP_START" -eq 1 ]; then
    export CLT_REPS=10000
  else
    export CLT_REPS=100000
  fi
  for worker in CLT_gaussian_hpc CLT_gaussian_small_n_hpc; do
    sbatch --array=1-100 --partition="$PARTITION" --account="$ACCOUNT" \
      --output="$CLT_RUN_ROOT/logs/%x-%A_%a.out" --error="$CLT_RUN_ROOT/logs/%x-%A_%a.err" \
      "simulation/ET_round2/${worker}.slurm"
  done
done
```

Once all four arrays complete, run the combined summary command above.
The following scheduler-free test uses the same two R worker entry points.
It computes all 196 bump integrals for four evaluation fits per size and
four independent reference fits at $n=500$, then writes `clt_gaussian_all_projections_R4_*.pdf`
and `jarque_bera_first4.csv`.

```bash
trial="$(mktemp -d)"
Rscript --vanilla simulation/ET_round2/CLT_gaussian_hpc.R \
  --mode=shard --reps=4 --num-shards=1 --shard-id=1 --shard-dir="$trial/shards"
Rscript --vanilla simulation/ET_round2/CLT_gaussian_small_n_hpc.R \
  --mode=shard --reps=4 --num-shards=1 --shard-id=1 --shard-dir="$trial/small-n-shards"
Rscript --vanilla simulation/ET_round2/CLT_gaussian_small_n_summary.R \
  --reps=4 --base-reps=4 --base-num-shards=1 --mc-block-size=4 \
  --base-shard-dir="$trial/shards" --small-n-shard-dir="$trial/small-n-shards" \
  --output-dir="$trial/result" --plot-dir="$trial/plot"
```

To run production without Slurm, use the same two worker commands with
`--num-shards=100` for each `--shard-id` from 1 through 100 in each
batch. Use `--reps=10000 --rep-start=1` for the base and
`--reps=100000 --rep-start=10001` for the continuation. Point both batches to the production
shard directories used by the combined summary.

## Original 10,000-replication run

The original target was 10,000 evaluation replications at each
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
export REPO_ROOT=/hpc/home/yx306/OT/UOT
export R_BIN=/opt/apps/rhel9/R-4.4.3/bin
export VENV_PYTHON=/hpc/home/yx306/venvs/uot/bin/python
export CLT_RUN_ROOT=/cwork/yx306/UOT/out/clt-gaussian-10000
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
your cluster requires a different request. In this configuration all
production logs, shards, merged RDS/CSV outputs, and PDF plots are under
`/cwork/yx306/UOT/out`. That directory must exist before `sbatch`, because
SLURM opens its log files before the script starts.

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

To summarize a completed production run directly from its immutable shards,
without creating the canonical merged checkpoints, run:

```bash
Rscript --vanilla simulation/ET_round2/CLT_gaussian_summary.R \
  --reps=10000 --rep-start=1 --num-shards=100 \
  --shard-dir="$CLT_SHARD_DIR" --output-dir="$CLT_RESULT_DIR" \
  --plot-dir="$CLT_PLOT_DIR"
```

The summary driver validates the complete design and disjoint seed streams
before writing its report, tables, covariance summary, and seven PDF figures.
Supply writable output and plot directories with the flags above.

### Small-sample extension

To add 10,000 replications only at `n=10,20`, preflight and submit the
evaluation-only array:

```bash
export CLT_SMALL_N_SHARD_DIR="$CLT_RUN_ROOT/small-n-shards"

Rscript --vanilla simulation/ET_round2/CLT_gaussian_small_n_hpc.R \
  --mode=preflight --shard-dir="$CLT_SMALL_N_SHARD_DIR"

sbatch --array=1-100 simulation/ET_round2/CLT_gaussian_small_n_hpc.slurm
```

After all 100 tasks finish, validate and combine the new sample sizes with the
existing run:

```bash
Rscript --vanilla simulation/ET_round2/CLT_gaussian_small_n_summary.R \
  --reps=10000 --base-reps=10000 --base-num-shards=100 \
  --base-shard-dir="$CLT_SHARD_DIR" \
  --small-n-shard-dir="$CLT_SMALL_N_SHARD_DIR" \
  --output-dir="$CLT_RESULT_DIR" --plot-dir="$CLT_PLOT_DIR"
```

The combined design is `n={10,20,50,100,200,500}` with 10,000 evaluation
replications per sample size and the original independent 10,000-run reference
at `n=500`. The uncertainty band uses 50 independent null blocks of size 400,
scaled to the 10,000-replication covariance comparison. It is an empirical
Monte Carlo diagnostic rather than a formal confidence interval.

To extend the complete six-size design from 10,000 to 100,000 replications,
run only IDs 10,001--100,000 in both worker arrays:

```bash
export CLT_REPS=100000
export CLT_REP_START=10001
export CLT_NUM_SHARDS=100

sbatch --array=1-100 simulation/ET_round2/CLT_gaussian_hpc.slurm
sbatch --array=1-100 simulation/ET_round2/CLT_gaussian_small_n_hpc.slurm
```

After both arrays complete, combine the base and continuation shard sets:

```bash
Rscript --vanilla simulation/ET_round2/CLT_gaussian_small_n_summary.R \
  --reps=100000 --base-reps=10000 \
  --base-num-shards=100 --extension-num-shards=100 \
  --extension-rep-start=10001 \
  --base-shard-dir="$CLT_SHARD_DIR" \
  --small-n-shard-dir="$CLT_SMALL_N_SHARD_DIR" \
  --output-dir="$CLT_RESULT_DIR" --plot-dir="$CLT_PLOT_DIR" \
  --mc-block-size=2000
```

The 100,000-run uncertainty calibration uses 100 independent `n=500` null
blocks of size 2,000. Base and continuation shards are validated separately
before their exact replication-ID ranges are combined.

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

For the original four-size, 10,000-replication run, the merged evaluation checkpoint must contain
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
