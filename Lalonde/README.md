# LaLonde earnings analysis

This example reproduces the UOT entry in Table 1 and the two LaLonde sensitivity plots. It uses 297 treated NSW participants, 425 randomized NSW controls for the experimental benchmark, and 2,871 observational PSID control rows. The PSID sample stacks the 2,490 base rows with its subsets of 253 and 128, as in the manuscript. The analysis downloads the five source files from the [Dehejia–Wahba data collection](https://users.nber.org/~rdehejia/nswdata2.html), checks their SHA-256 hashes, and caches them in `Lalonde/data/`.

## Run

From the repository root, create a Python environment with NumPy and SciPy. The versions in the requirements file were used with Python 3.12. The plots use R and ggplot2 (4.0 or later).

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r Lalonde/requirements.txt
export UOT_PYTHON="$PWD/.venv/bin/python"
Rscript -e 'install.packages(c("ggplot2", "MatchIt", "sbw", "quadprog", "MASS", "nnet", "numDeriv", "glmnet"), repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/CBPS/CBPS_0.23.tar.gz", repos=NULL, type="source")'
Rscript --vanilla Lalonde/Lalonde.R
```

The full Table 1 analysis uses archived CBPS 0.23. With the public ATT specification, version 0.23 gives −$404.72 and version 0.24 gives −$465.30. The full run was checked with R 4.6.1, MatchIt 4.8.2, sbw 1.2, quadprog 1.5.8, and ggplot2 4.0.3.

The full Table 1 command recomputes the grid and both figures, fits the comparison estimators, and writes Table 1. To compute only the UOT, experimental, and unadjusted entries, run `Rscript --vanilla Lalonde/Lalonde.R --uot-only`. To regenerate just the sensitivity grid and plots, run `Rscript --vanilla Lalonde/UOT_sensi.R`. The grid defaults to two worker processes. For example, `Rscript --vanilla Lalonde/UOT_sensi.R --jobs 3` uses three. Both R entry points use `python3` on the search path unless `UOT_PYTHON` selects another interpreter.

## Results

- [Table 1](results/table1.csv) uses the grid point that minimizes the absolute difference between the UOT-weighted ATT estimate and the experimental ATT estimate. The selected estimate is **$872.75** at ε = 0.000268 and ρ = 0.0001. The benchmark is **$886.30**. This is a retrospective selection using the experimental outcome.
- [Figure 4](results/UOT_matched_pairs_sensi.png) shows matching concentration. For each treated participant, controls are ranked by pair weight. The figure reports the median number needed to carry at least 95% of a participant's transported mass. The median is **2** at the selected setting and ranges from **1 to 1,434** on the grid.
- [Figure 5](results/UOT_ATT_sensi.png) shows the UOT-weighted ATT estimates for 1978 earnings. The UOT estimate weights treated participants by their transported row masses. The experimental ATT estimate weights them equally. The red box marks the Table 1 selection.
- The [grid CSV](results/converged_grid.csv) contains all 300 settings, weighted group means, estimates, transported mass, median 95% counts, optimality residuals, and the selected-row flag.

The solver uses seven pooled standardized covariates, Euclidean distance divided by the largest cross-group distance, and uniform input weights within each group. It solves the KL transport problem using log-domain Sinkhorn updates or a damped Newton solve for the treated potentials, depending on the ratio of the penalties. For each cell it computes both marginal masses in log space and checks the first-order conditions of the full primal objective. It stops if the maximum absolute residual exceeds 1e-10. Earnings enter only after the weights are fitted.
