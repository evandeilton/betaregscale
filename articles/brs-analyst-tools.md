# Analyst Tools for betaregscale

## Overview

This vignette presents analyst-oriented tools added to **betaregscale**
for post-estimation workflows:

- [`brs_table()`](https://evandeilton.github.io/betaregscale/reference/brs_table.md):
  model comparison tables.
- [`brs_marginaleffects()`](https://evandeilton.github.io/betaregscale/reference/brs_marginaleffects.md):
  average marginal effects (AME).
- [`brs_predict_scoreprob()`](https://evandeilton.github.io/betaregscale/reference/brs_predict_scoreprob.md):
  probabilities on the original score scale.
- [`autoplot.brs()`](https://evandeilton.github.io/betaregscale/reference/autoplot.brs.md):
  `ggplot2` diagnostics.
- [`brs_cv()`](https://evandeilton.github.io/betaregscale/reference/brs_cv.md):
  repeated k-fold cross-validation.

The examples use bounded scale data under the interval-censored beta
regression framework of Lopes (2024).

``` r
library(betaregscale)
```

## Mathematical foundations

### Complete likelihood by censoring type

For each observation $i$, let $\delta_{i} \in \{ 0,1,2,3\}$ indicate
exact, left-censored, right-censored, or interval-censored status. With
beta CDF $F( \cdot )$, beta density $f( \cdot )$, and interval endpoints
$\left\lbrack l_{i},u_{i} \right\rbrack$, the contribution is:

$$L_{i}(\theta) = \begin{cases}
{f\left( y_{i};a_{i},b_{i} \right),} & {\delta_{i} = 0,} \\
{F\left( u_{i};a_{i},b_{i} \right),} & {\delta_{i} = 1,} \\
{1 - F\left( l_{i};a_{i},b_{i} \right),} & {\delta_{i} = 2,} \\
{F\left( u_{i};a_{i},b_{i} \right) - F\left( l_{i};a_{i},b_{i} \right),} & {\delta_{i} = 3.}
\end{cases}$$

This is the basis for fitting, prediction, and validation metrics.

### Model-comparison metrics

[`brs_table()`](https://evandeilton.github.io/betaregscale/reference/brs_table.md)
reports:

- $\log L\left( \widehat{\theta} \right)$,
- $AIC = - 2\log L\left( \widehat{\theta} \right) + 2k$,
- $BIC = - 2\log L\left( \widehat{\theta} \right) + k\log n$,

where $k$ is the number of estimated parameters and $n$ is sample size.

### Average marginal effects (AME)

[`brs_marginaleffects()`](https://evandeilton.github.io/betaregscale/reference/brs_marginaleffects.md)
computes AME by finite differences:

$${AME}_{j} = \frac{1}{n}\sum\limits_{i = 1}^{n}\frac{{\widehat{g}}_{i}\left( x_{ij} + h \right) - {\widehat{g}}_{i}\left( x_{ij} \right)}{h},$$

with ${\widehat{g}}_{i}$ on the requested prediction scale (`response`
or `link`). For binary covariates $x_{j} \in \{ 0,1\}$, it uses the
discrete contrast
$\widehat{g}\left( x_{j} = 1 \right) - \widehat{g}\left( x_{j} = 0 \right)$.

### Score-scale probabilities

For integer scores $s \in \{ 0,\ldots,K\}$,
[`brs_predict_scoreprob()`](https://evandeilton.github.io/betaregscale/reference/brs_predict_scoreprob.md)
computes:

$$P(Y = s) = \begin{cases}
{F\left( {lim}/K \right),} & {s = 0,} \\
{1 - F\left( \left( K - {lim} \right)/K \right),} & {s = K,} \\
{F\left( \left( s + {lim} \right)/K \right) - F\left( \left( s - {lim} \right)/K \right),} & {1 \leq s \leq K - 1.}
\end{cases}$$

These probabilities are directly aligned with interval geometry on the
original instrument scale.

### Cross-validation log score

In
[`brs_cv()`](https://evandeilton.github.io/betaregscale/reference/brs_cv.md),
fold-level predictive quality includes:

$${log\_ score} = \frac{1}{n_{test}}\sum\limits_{i \in test}\log\left( p_{i} \right),$$

where $p_{i}$ is the predictive contribution from the same
censoring-rule piecewise definition shown above.

It also reports:

$${RMSE}_{yt} = \sqrt{\frac{1}{n_{test}}\sum\left( y_{t,i} - {\widehat{\mu}}_{i} \right)^{2}},\qquad{MAE}_{yt} = \frac{1}{n_{test}}\sum\left| y_{t,i} - {\widehat{\mu}}_{i} \right|.$$

## Reproducible workflow

### 1) Simulate data and fit candidate models

``` r
set.seed(2026)
n <- 220
d <- data.frame(
  x1 = rnorm(n),
  x2 = rnorm(n),
  z1 = rnorm(n)
)

sim <- brs_sim(
  formula = ~ x1 + x2 | z1,
  data = d,
  beta = c(0.20, -0.45, 0.25),
  zeta = c(0.15, -0.30),
  ncuts = 100,
  repar = 2
)

fit_fixed <- brs(y ~ x1 + x2, data = sim, repar = 2)
fit_var <- brs(y ~ x1 + x2 | z1, data = sim, repar = 2)
```

### 2) Compare models in one table

``` r
tab <- brs_table(
  fixed = fit_fixed,
  variable = fit_var,
  sort_by = "AIC"
)
tab
#>      model nobs npar    logLik      AIC      BIC pseudo_r2 exact left right
#> 1 variable  220    5 -931.8646 1873.729 1890.697    0.1381     0    8    22
#> 2    fixed  220    4 -939.3243 1886.649 1900.223    0.1381     0    8    22
#>   interval prop_exact prop_left prop_right prop_interval
#> 1      190          0    0.0364        0.1        0.8636
#> 2      190          0    0.0364        0.1        0.8636
```

### 3) Estimate average marginal effects

``` r
me_mean <- brs_marginaleffects(
  fit_var,
  model = "mean",
  type = "response",
  interval = TRUE,
  n_sim = 120
)
me_mean
#>   variable         ame  std.error    ci.lower    ci.upper model     type   n
#> 1       x1 -0.10350449 0.01859635 -0.13287991 -0.06613533  mean response 220
#> 2       x2  0.07342329 0.01729846  0.04748093  0.10839857  mean response 220

me_precision <- brs_marginaleffects(
  fit_var,
  model = "precision",
  type = "link",
  interval = FALSE
)
me_precision
#>   variable        ame std.error ci.lower ci.upper     model type   n
#> 1       z1 -0.3360715        NA       NA       NA precision link 220
```

### 4) Predict score probabilities

``` r
P <- brs_predict_scoreprob(fit_var, scores = 0:10)
dim(P)
#> [1] 220  11
head(P[, 1:5])
#>          score_0     score_1     score_2     score_3     score_4
#> [1,] 0.014623303 0.017790410 0.014558715 0.013036068 0.012076830
#> [2,] 0.021451858 0.022805666 0.017769488 0.015484433 0.014073831
#> [3,] 0.046547692 0.031830606 0.021590702 0.017455271 0.015064658
#> [4,] 0.073078453 0.037104253 0.023313719 0.018090062 0.015175980
#> [5,] 0.007646182 0.010508656 0.009021869 0.008303000 0.007845942
#> [6,] 0.005023701 0.005045276 0.003869686 0.003351163 0.003037766
summary(rowSums(P))
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> 0.009836 0.077558 0.133484 0.150893 0.212978 0.526155
```

`rowSums(P)` should be close to 1 (up to numerical tolerance and score
subset selection).

### 5) ggplot2 diagnostics

``` r
autoplot.brs(fit_var, type = "calibration")
```

![](brs-analyst-tools_files/figure-html/autoplot-1.png)

``` r
autoplot.brs(fit_var, type = "score_dist", scores = 0:20)
```

![](brs-analyst-tools_files/figure-html/autoplot-2.png)

``` r
autoplot.brs(fit_var, type = "cdf", max_curves = 4)
```

![](brs-analyst-tools_files/figure-html/autoplot-3.png)

``` r
autoplot.brs(fit_var, type = "residuals_by_delta", residual_type = "rqr")
```

![](brs-analyst-tools_files/figure-html/autoplot-4.png)

### 6) Repeated k-fold cross-validation

``` r
cv_res <- brs_cv(
  y ~ x1 + x2 | z1,
  data = sim,
  k = 3,
  repeats = 1,
  repar = 2,
  seed = 2026
)

cv_res
#>   repeat fold n_train n_test log_score   rmse_yt    mae_yt converged error
#> 1      1    1     146     74 -4.307714 0.3214558 0.2861862      TRUE  <NA>
#> 2      1    2     147     73 -4.317109 0.3448886 0.3031531      TRUE  <NA>
#> 3      1    3     147     73 -4.167553 0.3454093 0.3021977      TRUE  <NA>
colMeans(cv_res[, c("log_score", "rmse_yt", "mae_yt")], na.rm = TRUE)
#>  log_score    rmse_yt     mae_yt 
#> -4.2641252  0.3372512  0.2971790
```

## Practical interpretation

- Prefer the model with lower AIC/BIC and better predictive `log_score`.
- Use AME on the response scale to communicate expected change in mean
  score (on the unit interval) from small covariate shifts.
- Use score probabilities to translate model outputs back to clinically
  interpretable scale categories.
- Inspect calibration and residual-by-censoring plots before inference.

## References

- Lopes, J. E. (2024). *Beta Regression for Interval-Censored
  Scale-Derived Outcomes*. MSc Dissertation, PPGMNE/UFPR.
- Ferrari, S. and Cribari-Neto, F. (2004). Beta regression for modelling
  rates and proportions. *Journal of Applied Statistics*, 31(7),
  799-815.
- Hawker, G. A., Mian, S., Kendzerska, T., and French, M. (2011).
  Measures of adult pain: VAS, NRS, MPQ, SF-MPQ, CPGS, SF-36 BPS, and
  ICOAP. *Arthritis Care and Research*, 63(S11), S240-S252.
- Hjermstad, M. J., Fayers, P. M., Haugen, D. F., et al. (2011).
  Comparing NRS, VRS, and VAS pain scales in adults: a systematic
  review. *Journal of Pain and Symptom Management*, 41(6), 1073-1093.
