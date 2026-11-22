# UOT

Code for simulation studies and empirical examples around unbalanced optimal transport (UOT), including synthetic experiments, Lalonde-based treatment effect analyses, and ACS-based spouse-matching experiments.

## Repository layout

- `simulation/`
  [Figure 1](simulation/ET_round2/figure1/) compares the normalized empirical and population Gaussian couplings using the same squared Euclidean cost. From the repository root, run `Rscript --vanilla simulation/ET_round2/MMD_heatmap.R` with MASS and ggplot2 to regenerate the 2,450 fits and both panels.
  [Figure 2](simulation/ET_round2/figure2/) shows raw matched mass beyond Euclidean distance 0.1 as both KL penalties decrease, alongside the bound from the proof of Proposition 3. Run `Rscript --vanilla simulation/MR_R1Q2_2.R` from the repository root to reproduce the figure and the 300 fit diagnostics.
  `Major_revision.R` contains the Gaussian UOT examples and moment comparisons.
  `MR_*.R`, `One_Rep_*.R`, and `Subsampling_method.R` contain additional simulation workflows.
  `simulation/ET_round2/` contains the revised Gaussian UOT CLT experiment,
  including SLURM instructions for the six-sample-size study with up to
  100,000 replications per sample size and an independent reference batch.

- [`Lalonde/`](Lalonde/README.md)
  Replicates Table 1 and Figures 4–5 with convergence-checked KL transport weights. `Lalonde.R` runs the full example and `UOT_sensi.R` regenerates the tuning grid and plots.

- [`ACS/`](ACS/README.md)
  Reproduces Figure 6 from paired ages in the 2024 ACS extract. The driver checks convergence of each product-reference KL fit and saves the cell summaries and fit diagnostics.

- `R/`
  Shared R utilities for transport calculations, plotting, matching summaries, and analytic examples.
  `thresh_ot_func.R` is the main helper file used across scripts.

- `python/`
  Python code called from R through `reticulate`.
  `sinkhorn_unbalanced_tv.py` exposes the unbalanced Sinkhorn routines used by the R scripts.

- `ATO.R` and `analyze_allATO.R`
  Older simulation code from the earlier arXiv workflow. These scripts are still useful if you want the previous ATO-style Monte Carlo setup.

## Requirements

For the LaLonde example, follow its [setup and run instructions](Lalonde/README.md). The remaining R/Python workflows use the following tools.

- R with packages used by the scripts, including `reticulate`, `ggplot2`, `dplyr`, `cobalt`, `sandwich`, and `lmtest`
- Python accessible from `reticulate`
- The Python OT package `POT` available in that Python environment

Several scripts currently hardcode a Python path such as:

```r
use_python("~/Library/r-miniconda-arm64/bin/python", required = TRUE)
```

Change that line to match your local Python environment before running the code.

## Data inputs

- The LaLonde entry points download and verify the five NSW and PSID files from NBER on first use. Subsequent runs use the cached files in `Lalonde/data/`.
- The Figure 6 driver reads `data/usa_00002.csv.gz` or an extract path supplied on the command line.
- Some simulation scripts expect local saved objects such as `rng_state.rds` or output folders like `result/` and `plot/`. Create those if you want to reproduce those exact runs.

## Typical entry points

Run these from the repository root.

```bash
Rscript Lalonde/Lalonde.R
Rscript Lalonde/UOT_sensi.R
Rscript --vanilla ACS/figure6.R
Rscript simulation/Major_revision.R
```

For the revised CLT experiment and its HPC workflow, see
[`simulation/ET_round2/README.md`](simulation/ET_round2/README.md). Start with
the environment preflight:

```bash
Rscript --vanilla simulation/ET_round2/CLT_gaussian_hpc.R --mode=preflight
```

For the older simulation pipeline:

```bash
Rscript analyze_allATO.R
```

Edit the top-level parameters inside each script before long runs. Most scripts are written as research workflows rather than packaged functions, so sample sizes, tuning parameters, file paths, and output locations are configured inline.

## Notes

- This repository is organized as a research codebase, not an R package.
- Many scripts generate figures or intermediate objects interactively and may assume folders already exist.
- If you want a fully reproducible run on a new machine, start by fixing the Python path, installing the required R/Python dependencies, and checking that the expected local data files and output directories exist.
