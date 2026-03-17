# UOT

Code for simulation studies and empirical examples around unbalanced optimal transport (UOT), including synthetic experiments, Lalonde-based treatment effect analyses, and ACS-based spouse-matching experiments.

## Repository layout

- `simulation/`
  Main synthetic experiments and revision scripts.
  `Major_revision.R` contains the Gaussian UOT examples and moment comparisons.
  `MR_*.R`, `One_Rep_*.R`, and `Subsampling_method.R` contain additional simulation workflows.

- `Lalonde/`
  Empirical treatment-effect scripts based on the Lalonde / NSW / PSID data.
  `Lalonde.R` runs the main UOT analysis and benchmark comparisons.
  `UOT_sensi.R` performs a grid search over UOT tuning parameters.
  `lalonde_help.R` contains helper functions for loading data and evaluating matches.

- `ACS/`
  ACS preprocessing and bootstrap experiments for spouse-pair matching.
  `ACS1_preprocess_data.R` reads ACS microdata, builds spouse pairs, and runs one bootstrap/UOT pipeline.
  `ACS2_bootstrap.R` contains the bootstrap workflow used after preprocessing.
  `ACS_help.R` provides helper utilities for ACS data cleaning and pair construction.

- `R/`
  Shared R utilities for transport calculations, plotting, matching summaries, and analytic examples.
  `thresh_ot_func.R` is the main helper file used across scripts.

- `python/`
  Python code called from R through `reticulate`.
  `sinkhorn_unbalanced_tv.py` exposes the unbalanced Sinkhorn routines used by the R scripts.

- `ATO.R` and `analyze_allATO.R`
  Older simulation code from the earlier arXiv workflow. These scripts are still useful if you want the previous ATO-style Monte Carlo setup.

## Requirements

- R with packages used by the scripts, including `reticulate`, `ggplot2`, `dplyr`, `cobalt`, `sandwich`, and `lmtest`
- Python accessible from `reticulate`
- The Python OT package `POT` available in that Python environment

Several scripts currently hardcode a Python path such as:

```r
use_python("~/Library/r-miniconda-arm64/bin/python", required = TRUE)
```

Change that line to match your local Python environment before running the code.

## Data inputs

- `Lalonde/Lalonde.R` and `Lalonde/UOT_sensi.R` download the NSW and PSID files directly from the Dehejia-Wahba / NBER source.
- `ACS/ACS1_preprocess_data.R` expects ACS microdata under `./data/`, currently reading `./data/usa_00002.csv.gz`.
- Some simulation scripts expect local saved objects such as `rng_state.rds` or output folders like `result/` and `plot/`. Create those if you want to reproduce those exact runs.

## Typical entry points

Run these from the repository root.

```bash
Rscript Lalonde/Lalonde.R
Rscript Lalonde/UOT_sensi.R
Rscript ACS/ACS1_preprocess_data.R
Rscript simulation/Major_revision.R
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
