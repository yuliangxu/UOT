# Figure 6 · ACS spouse age matching

From the repository root, run the following with R, ggplot2 (4.0.3 used for the manuscript) and scales (1.4.0). The solver uses base R.

```sh
Rscript --vanilla ACS/figure6.R data/usa_00002.csv.gz ACS/results
```

The inputs are the 2024 IPUMS USA ACS extract and its household spouse links. The driver reads YEAR, SERIAL, PERNUM, SPLOC, SEX and AGE. The manuscript extract (SHA-256 `37fdb3dd5fd8a8357dea950adf307283a9db8376a0c364027a4c77d66130bf92`) yields 774,014 linked male–female couples. Figure 6 uses all of them. Pooled ages span 15–96 years. Sampled couple IDs are ordered as in the original bootstrap, and each draw contributes one age to each histogram. The analysis uses 30 equal-width bins and unweighted draws of 1,000 couples, with seeds 124–173 reused at every parameter setting.

An earlier version excluded one couple because the husband's recorded age of 15 raised concern about a possible reporting error. The extract has 21 husbands recorded as 15, including two whose wives are recorded as 20. The recorded ages gave no basis for excluding just one of those two couples. The current analysis applies the same spouse-link criteria to every pair.

The distance between bin centers is divided by its maximum. The solver uses the product of the two empirical histograms as the KL reference, restricts its calculation to occupied bins, and restores zero entries on the full grid. It checks the full first-order conditions in log space and the total-mass bound before returning each fit. A failed check stops the analysis.

- `acs.png` shows the mean and the pointwise 2.5% and 97.5% quantiles across 50 bootstrap plans at each of 12 parameter settings. The shared power-transformed color scale is labeled in transport-mass units.
- `figure6_panels.csv` contains the values plotted in all 36 panels. Male bins index rows and female bins index columns, both in increasing age order.
- `figure6_histograms.csv` records the paired bootstrap marginals on the original age breaks. `figure6_diagnostics.csv` records the occupied-bin counts, solver iterations, maximum absolute first-order residual and total mass for all 600 fits.

In the mean panel at ε = 0.1 and ρ = 0.01, the total matched mass is 0.30460389 and the largest cell has mass 0.00513165. All 600 fits pass the 1e-10 first-order residual tolerance. An independent solver agreed on all 32,400 plotted mean and quantile values to within 1.1e-12.
