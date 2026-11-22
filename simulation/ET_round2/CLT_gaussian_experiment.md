# Numerical illustration of the Gaussian-process limit

## Experiment setting

We illustrate the asymptotic distribution in Theorem 1 using two Gaussian
marginals in dimension $d=2$:

$$
\mu_0=N_2(a,A), \qquad \mu_1=N_2(b,B),
$$

where

$$
a=(-0.4,0.1)^\top, \qquad b=(0.5,-0.3)^\top,
$$

and

$$
A=
\begin{pmatrix}
0.7 & -0.5\\
-0.5 & 1
\end{pmatrix},
\qquad
B=
\begin{pmatrix}
1 & -0.34\\
-0.34 & 0.7
\end{pmatrix}.
$$

The UOT problem uses the fixed squared Euclidean cost
$c(x^{(0)},x^{(1)})=\lVert x^{(0)}-x^{(1)}\rVert^2$, entropic parameter
$\varepsilon=0.5$, and marginal relaxation parameter $\rho=1$. The empirical
solver uses KL regularization relative to the product of the empirical
marginals. The transport plan is not normalized after optimization.

For these parameters, the population UOT coupling has the closed form

$$
\gamma_0=m_\gamma N_4(\mu_\gamma,H_\gamma),
$$

with transported mass $m_\gamma=0.588149$ and mean

$$
\mu_\gamma=
(-0.131862,-0.118713,0.132774,-0.150429)^\top.
$$

For each $n\in\{10,20,50,100,200,500\}$, we independently draw

$$
X_i^{(0)}\sim\mu_0,\quad i=1,\ldots,n,
\qquad
X_j^{(1)}\sim\mu_1,\quad j=1,\ldots,n.
$$

Here, the superscripts $(0)$ and $(1)$ identify the two marginals, while the
subscripts $i$ and $j$ identify observations. In particular,
$X_i^{(0)}=(X_{i,1}^{(0)},X_{i,2}^{(0)})$ and
$X_j^{(1)}=(X_{j,1}^{(1)},X_{j,2}^{(1)})$ are observed random vectors, not
population parameters. Each empirical coupling support point is

$$
Z_{ij}=\bigl(X_i^{(0)},X_j^{(1)}\bigr)\in\mathbb R^4.
$$

We compute the empirical UOT coupling $\widehat\gamma_n$ on these support
points. The replication counts for the original diagnostic suite and the
high-precision covariance study are given below. For any fixed test function
$f$, define its empirical and population coupling integrals by

$$
\widehat\theta_n(f)=\int f\,d\widehat\gamma_n,
\qquad
\theta_0(f)=\int f\,d\gamma_0.
$$

The **scaled fluctuation** associated with $f$ is

$$
\Delta_n(f)
=
\sqrt n\left\{\widehat\theta_n(f)-\theta_0(f)\right\}
=
\sqrt n
\left\{
\int f\,d\widehat\gamma_n-
\int f\,d\gamma_0
\right\}.
$$

Here, “scaled” refers only to multiplication by the central-limit rate
$\sqrt n$; it does not mean division by an estimated standard deviation.

We use the common per-marginal sample size $n$ as the asymptotic index. If the
theorem's $N$ is instead instantiated as the total sample size
$N_{\mathrm{tot}}=2n$ in this equal-sample design, its scaled statistic is
$\sqrt{2}\,\Delta_n(f)$. This constant convention changes the limiting
covariance by a factor of two, but it does not change Gaussianity,
correlations, standardized Q--Q plots, or relative covariance comparisons.

Raw transport-matrix entries cannot be compared across replications because
their support points change with each sample. To make the test functions
distinct from the observed samples, write a generic point in the population
coupling as

$$
z=\bigl(z^{(0)},z^{(1)}\bigr),
\qquad
z^{(0)}=(z_1^{(0)},z_2^{(0)}),
\qquad
z^{(1)}=(z_1^{(1)},z_2^{(1)}).
$$

For every coordinate pair $(r,s)\in\{1,2\}^2$, we use a fixed $7\times7$
collection of bounded Gaussian bump functions on the projection
$(z_r^{(0)},z_s^{(1)})$. For fixed grid centers $u$ and $v$, define

$$
f_{u,v}^{(r,s)}(z)=
\exp\left[
-\frac{(z_r^{(0)}-u)^2+(z_s^{(1)}-v)^2}{2h_{r,s}^2}
\right].
$$

The four families cover the two matched-coordinate projections and both
cross-coordinate projections:

| Projection | Type | Bandwidth |
|---|---|---:|
| $(z_1^{(0)},z_1^{(1)})$ | matched $(1,1)$ | $h_{1,1}=0.4157$ |
| $(z_1^{(0)},z_2^{(1)})$ | cross $(1,2)$ | $h_{1,2}=0.3973$ |
| $(z_2^{(0)},z_1^{(1)})$ | cross $(2,1)$ | $h_{2,1}=0.4381$ |
| $(z_2^{(0)},z_2^{(1)})$ | matched $(2,2)$ | $h_{2,2}=0.4197$ |

The centers and bandwidths are fixed from the closed-form population coupling
before any empirical samples are drawn. Each axis has seven equally spaced
centers spanning its coordinate mean under the normalized population coupling
$\gamma_0/m_\gamma$, plus or minus two corresponding standard deviations. If
$\sigma_r^{(0)}$ and $\sigma_s^{(1)}$ denote these projected population
standard deviations, then

$$
h_{r,s}=\frac{1}{2}\left(\frac{\sigma_r^{(0)}+\sigma_s^{(1)}}{2}\right)
=\frac{\sigma_r^{(0)}+\sigma_s^{(1)}}{4}.
$$

For each family, $u$ ranges over the source-coordinate grid $\mathcal U_r$
and $v$ over the target-coordinate grid $\mathcal V_s$. The four grids are

$$
\begin{aligned}
\mathcal U_1={}&\{-1.702,-1.179,-0.655,-0.132,0.391,0.915,1.438\},\\
\mathcal U_2={}&\{-1.868,-1.285,-0.702,-0.119,0.464,1.047,1.631\},\\
\mathcal V_1={}&\{-1.623,-1.038,-0.453,0.133,0.718,1.303,1.889\},\\
\mathcal V_2={}&\{-1.759,-1.223,-0.687,-0.150,0.386,0.922,1.458\}.
\end{aligned}
$$

Thus, the experiment contains $4\times49=196$ bump functions. For projection
$(r,s)$, the scaled local fluctuation at grid location $(u,v)$ is

$$
\Delta_n^{(r,s)}(u,v)
=
\Delta_n\left(f_{u,v}^{(r,s)}\right)
=
\sqrt n
\left\{
\int f_{u,v}^{(r,s)}\,d\widehat\gamma_n-
\int f_{u,v}^{(r,s)}\,d\gamma_0
\right\}.
$$

For a fixed projection, the 49 functions differ only through their bump
centers $(u,v)$. The **central bump** is the member whose center is the middle
grid point,

$$
u=\frac{1}{m_\gamma}\int z_r^{(0)}\,d\gamma_0=(\mu_\gamma)_r,
\qquad
v=\frac{1}{m_\gamma}\int z_s^{(1)}\,d\gamma_0=(\mu_\gamma)_{2+s}.
$$

The four central-bump locations are

| Projection | Central location $(u,v)$ |
|---|---:|
| $(1,1)$ | $(-0.131862,\phantom{-}0.132774)$ |
| $(1,2)$ | $(-0.131862,-0.150429)$ |
| $(2,1)$ | $(-0.118713,\phantom{-}0.132774)$ |
| $(2,2)$ | $(-0.118713,-0.150429)$ |

These are coordinate centers of the population UOT coupling, not generally
the means $a$ and $b$ of the two input marginals. Here, “central” describes
the location of the bump center; it does not mean that the function is
evaluated only at that point. In every replication, each
$f_{u,v}^{(r,s)}$ is evaluated at all empirical support points $Z_{ij}$ and
then integrated using all entries of the empirical transport plan. The other
48 functions in the same projection family are centered at noncentral
locations in that two-dimensional projected coordinate space.

We call this a local or smoothed projected-mass fluctuation because
$f_{u,v}^{(r,s)}$ emphasizes coupling mass near
$(z_r^{(0)},z_s^{(1)})=(u,v)$ while integrating over the other two
coordinates. The bumps overlap, so the tiles are not disjoint histogram bins
or pointwise coupling-density estimates. Each heatmap tile represents a
different test function, not a different argument of one fixed function.
These four projection families cover every source-target coordinate pair,
although they are still two-dimensional projections rather than arbitrary
four-dimensional test functions. Since $0<f_{u,v}^{(r,s)}\leq1$ and the
collection is finite and fixed, it is a uniformly bounded Donsker class of
test functions of the type allowed in Theorem 1.

Their population integrals are evaluated exactly under the Gaussian coupling.
We also include $f(z)=1$, whose integral is the total transported mass.

The theorem does not provide a closed-form covariance kernel for its limiting
Gaussian process. The experiment therefore has two layers. The original
diagnostic suite uses 2,000 evaluation replications at
$n\in\{50,100,200,500\}$ and a separate 2,000-run reference batch at $n=500$.
These runs produce the heatmaps, histograms, Q--Q plots, and reference
correlation plot described below.

The high-precision covariance study uses 100,000 evaluation replications at
each $n\in\{10,20,50,100,200,500\}$ and an independent 100,000-run reference
batch at $n=500$. Evaluation and reference batches use disjoint random-number
streams. The reference covariance $\widehat\Sigma_{\mathrm{ref}}$ remains a
finite-$n$ Monte Carlo proxy, not the analytic covariance kernel of the
limiting Gaussian process.

## Interpretation of the plots

Except for the high-precision covariance figure in the final subsection, the
plots below use the original 2,000-run design with
$n\in\{50,100,200,500\}$.

### Representative fluctuations

[`clt_gaussian_all_projections_representative_fluctuations.pdf`](plot/clt_gaussian_all_projections_representative_fluctuations.pdf)
shows one realization of the scaled fluctuation over all four projection grids
for each $n$. Red and blue regions represent positive and negative
fluctuations of the smoothed projected-mass functionals. Their magnitudes
remain nondegenerate after multiplication by $\sqrt{n}$, as expected under a
central limit theorem.

### Mean and standard-deviation heatmaps

[`clt_gaussian_all_projections_mean_fluctuations.pdf`](plot/clt_gaussian_all_projections_mean_fluctuations.pdf)
shows that the Monte Carlo mean of the scaled bump-function fluctuations is
generally smaller at larger $n$, with some Monte Carlo variation. Across the
196 plotted bumps, the maximum absolute mean is $0.0237$ at $n=50$ and
$0.0106$ at $n=500$. In the numerical summary that also includes total mass,
the corresponding maxima are $0.0689$ and $0.0160$.

[`clt_gaussian_all_projections_sd_fluctuations.pdf`](plot/clt_gaussian_all_projections_sd_fluctuations.pdf)
shows stable spatial variance patterns across sample sizes. The matched
projections exhibit diagonal structure, while the cross projections have
broader patterns reflecting dependence between different coordinates.
Stability after $\sqrt{n}$ scaling is consistent with a nondegenerate limiting
process.

### Histograms and Q--Q plots

[`clt_gaussian_all_projections_histograms.pdf`](plot/clt_gaussian_all_projections_histograms.pdf)
and
[`clt_gaussian_all_projections_qqplots.pdf`](plot/clt_gaussian_all_projections_qqplots.pdf)
examine five prespecified functionals:

- `total_mass` is the fluctuation in total transported mass.
- `p11_center` is the central bump for the matched $(1,1)$ projection.
- `p12_center` is the central bump for the cross $(1,2)$ projection.
- `p21_center` is the central bump for the cross $(2,1)$ projection.
- `p22_center` is the central bump for the matched $(2,2)$ projection.

Only these five functionals are shown as histograms to keep the figure
readable; the heatmaps contain all 196 bump-function fluctuations. Each Q--Q
panel therefore concerns one fixed function evaluated over repeated empirical
couplings. It does not combine the 49 bumps in a projection or evaluate one
function at 49 locations.

The red curves in the histograms are centered Gaussian densities whose
variances come only from the independent 2,000-run $n=500$ reference batch.
Overall, the histograms and Q--Q plots are broadly consistent with the
Gaussian reference. The clearest small-sample departure occurs for total
mass, which has a negative finite-sample bias and a heavier lower tail at
$n=50$; this discrepancy diminishes as $n$ grows.

### Covariance structure

[`clt_gaussian_all_projections_reference_correlation.pdf`](plot/clt_gaussian_all_projections_reference_correlation.pdf)
displays the correlation matrix estimated from the original 2,000-run
$n=500$ reference batch. The labeled blocks separate total mass and the four
projection families. Its structured positive and negative bands show
dependence both within and across projections, as expected for
finite-dimensional evaluations of one Gaussian process.

The high-precision figure
`clt_gaussian_n10_to_500_R100000_frobenius_distance_mc_band.pdf` compares each
100,000-run covariance estimate with the independent 100,000-run reference.
The relative Frobenius distance is

$$
D_n=
\frac{\lVert\widehat\Sigma_n-\widehat\Sigma_{\mathrm{ref}}\rVert_F}
{\lVert\widehat\Sigma_{\mathrm{ref}}\rVert_F}.
$$

| $n$ | $D_n$ | Monte Carlo classification |
|---:|---:|:---|
| 10 | 0.1134 | Above band |
| 20 | 0.0600 | Above band |
| 50 | 0.0291 | Above band |
| 100 | 0.0206 | Above band |
| 200 | 0.0172 | Within band |
| 500 | 0.0151 | Within band |

The empirical 95% Monte Carlo band is $[0.0133,0.0186]$, with RMS benchmark
$0.0158$. It is constructed from all pairwise covariance distances among 100
independent $n=500$ blocks of size 2,000, rescaled by

$$
\sqrt{\frac{2000-1}{100000-1}}.
$$

Because the pairwise block distances are dependent, this band is a Monte
Carlo diagnostic rather than a formal confidence interval. The monotone
decrease in $D_n$ shows a clear sample-size effect through $n=100$. At
$n=200$ and $n=500$, the remaining discrepancy is indistinguishable from
Monte Carlo covariance-estimation error at this resolution.

The same 100,000-run study gives mean transported masses $0.5475$, $0.5658$,
$0.5783$, $0.5830$, $0.5854$, and $0.5870$ for
$n=10,20,50,100,200,500$, respectively, approaching
$m_\gamma=0.5881$. Together, the centering, Gaussian diagnostics, and
high-precision covariance stabilization provide numerical support for the
finite-dimensional Gaussian limit in Theorem 1; they do not constitute a
formal joint-normality test of the full 197-dimensional vector.
