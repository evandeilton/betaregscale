# Introduction to betaregscale

## Overview

The **betaregscale** package provides maximum-likelihood estimation of
beta regression models for responses derived from bounded rating scales.
Common examples include pain intensity scales (NRS-11, NRS-21, NRS-101),
Likert-type scales, product quality ratings, and any instrument whose
response can be mapped to the open interval $(0,1)$.

The key idea is that a discrete score recorded on a bounded scale
carries measurement uncertainty inherent to the instrument. For
instance, a pain score of $y = 6$ on a 0–10 NRS is not an exact value
but rather represents a range: after rescaling to $(0,1)$, the
observation is treated as interval-censored in
$\lbrack 0.55,0.65\rbrack$. The package uses the beta distribution to
model such data, building a complete likelihood that supports mixed
censoring types within the same dataset.

## Installation

``` r
# Development version from GitHub:
# install.packages("remotes")
remotes::install_github("evandeilton/betaregscale")
```

``` r
library(betaregscale)
```

## Censoring types

The complete likelihood (Lopes, 2024, Eq. 2.24) supports four censoring
types, automatically classified by
[`check_response()`](https://evandeilton.github.io/betaregscale/reference/check_response.md):

| $\delta$ | Type                     | Likelihood contribution                                                       |
|:--------:|:-------------------------|:------------------------------------------------------------------------------|
|    0     | Exact (uncensored)       | $f\left( y_{i};\, a_{i},b_{i} \right)$                                        |
|    1     | Left-censored ($y = 0$)  | $F\left( u_{i};\, a_{i},b_{i} \right)$                                        |
|    2     | Right-censored ($y = K$) | $1 - F\left( l_{i};\, a_{i},b_{i} \right)$                                    |
|    3     | Interval-censored        | $F\left( u_{i};\, a_{i},b_{i} \right) - F\left( l_{i};\, a_{i},b_{i} \right)$ |

where $f( \cdot )$ and $F( \cdot )$ are the beta density and CDF,
$\left\lbrack l_{i},u_{i} \right\rbrack$ are the interval endpoints, and
$\left( a_{i},b_{i} \right)$ are the beta shape parameters derived from
$\mu_{i}$ and $\phi_{i}$ via the chosen reparameterization.

## Interval construction

Scale observations are mapped to $(0,1)$ with uncertainty intervals
controlled by the `type` argument:

- `"m"` (midpoint): $y_{t} = y/K$, interval
  $\left\lbrack y_{t} - 0.5/K,\; y_{t} + 0.5/K \right\rbrack$
- `"l"` (left-aligned): $y_{t} = y/K$, interval
  $\left\lbrack y_{t},\; y_{t} + 1/K \right\rbrack$
- `"r"` (right-aligned): $y_{t} = y/K$, interval
  $\left\lbrack y_{t} - 1/K,\; y_{t} \right\rbrack$

where $K$ is the number of scale categories (`ncuts`).

``` r
# Illustrate check_response with a 0-10 NRS scale
y_example <- c(0, 3, 5, 7, 10)
cr <- check_response(y_example, type = "m", ncuts = 10)
cr
#>         left   right      yt  y delta
#> [1,] 0.00001 0.05000 0.00001  0     1
#> [2,] 0.25000 0.35000 0.30000  3     3
#> [3,] 0.45000 0.55000 0.50000  5     3
#> [4,] 0.65000 0.75000 0.70000  7     3
#> [5,] 0.95000 0.99999 0.99999 10     2
```

The `delta` column shows that $y = 0$ is left-censored ($\delta = 1$),
$y = 10$ is right-censored ($\delta = 2$), and all interior values are
interval-censored ($\delta = 3$).

## Data preparation with `bs_prepare()`

In practice, analysts may want to supply their own censoring indicators
or interval endpoints rather than relying on the automatic
classification of
[`check_response()`](https://evandeilton.github.io/betaregscale/reference/check_response.md).
The
[`bs_prepare()`](https://evandeilton.github.io/betaregscale/reference/bs_prepare.md)
function provides a flexible, validated bridge between raw analyst data
and
[`betaregscale()`](https://evandeilton.github.io/betaregscale/reference/betaregscale.md).

It supports four input modes:

### Mode 1: Score only (automatic)

``` r
# Equivalent to check_response <U+2014> delta inferred from y
d1 <- data.frame(y = c(0, 3, 5, 7, 10), x1 = rnorm(5))
bs_prepare(d1, ncuts = 10)
#> bs_prepare: n = 5 | exact = 0, left = 1, right = 1, interval = 3
#>      left   right      yt  y delta           x1
#> 1 0.00001 0.05000 0.00001  0     1 -1.400043517
#> 2 0.25000 0.35000 0.30000  3     3  0.255317055
#> 3 0.45000 0.55000 0.50000  5     3 -2.437263611
#> 4 0.65000 0.75000 0.70000  7     3 -0.005571287
#> 5 0.95000 0.99999 0.99999 10     2  0.621552721
```

### Mode 2: Score + explicit censoring indicator

``` r
# Analyst specifies delta directly
d2 <- data.frame(
  y     = c(50, 0, 99, 50),
  delta = c(0, 1, 2, 3),
  x1    = rnorm(4)
)
bs_prepare(d2, ncuts = 100)
#> Warning: Observation(s) 3: delta = 2 (right-censored) but y != 100.
#> bs_prepare: n = 4 | exact = 1, left = 1, right = 1, interval = 1
#>      left   right      yt  y delta         x1
#> 1 0.50000 0.50000 0.50000 50     0  1.1484116
#> 2 0.00001 0.00500 0.00001  0     1 -1.8218177
#> 3 0.99500 0.99999 0.99999 99     2 -0.2473253
#> 4 0.49500 0.50500 0.50000 50     3 -0.2441996
```

### Mode 3: Interval endpoints with NA patterns

When the analyst provides `left` and/or `right` columns, censoring is
inferred from the NA pattern:

``` r
d3 <- data.frame(
  left  = c(NA, 20, 30, NA),
  right = c(5, NA, 45, NA),
  y     = c(NA, NA, NA, 50),
  x1    = rnorm(4)
)
bs_prepare(d3, ncuts = 100)
#> bs_prepare: n = 4 | exact = 1, left = 1, right = 1, interval = 1
#>    left   right    yt  y delta         x1
#> 1 1e-05 0.05000 0.025 NA     1 -0.2827054
#> 2 2e-01 0.99999 0.600 NA     2 -0.5536994
#> 3 3e-01 0.45000 0.375 NA     3  0.6289820
#> 4 5e-01 0.50000 0.500 50     0  2.0650249
```

### Mode 4: Analyst-supplied intervals

When the analyst provides `y`, `left`, and `right` simultaneously, their
endpoints are used directly (rescaled by $K$):

``` r
d4 <- data.frame(
  y     = c(50, 75),
  left  = c(48, 73),
  right = c(52, 77),
  x1    = rnorm(2)
)
bs_prepare(d4, ncuts = 100)
#> bs_prepare: n = 2 | exact = 0, left = 0, right = 0, interval = 2
#>   left right   yt  y delta         x1
#> 1 0.48  0.52 0.50 50     3 -1.6309894
#> 2 0.73  0.77 0.75 75     3  0.5124269
```

### Using prepared data with `betaregscale()`

Data processed by
[`bs_prepare()`](https://evandeilton.github.io/betaregscale/reference/bs_prepare.md)
is automatically detected by
[`betaregscale()`](https://evandeilton.github.io/betaregscale/reference/betaregscale.md)
— the internal
[`check_response()`](https://evandeilton.github.io/betaregscale/reference/check_response.md)
step is skipped:

``` r
set.seed(42)
n <- 200
dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
sim <- betaregscale_simulate(
  formula = ~ x1 + x2, data = dat,
  beta = c(0.2, -0.5, 0.3), phi = 1 / 5
)
prep <- bs_prepare(sim, ncuts = 100)
#> bs_prepare: n = 200 | exact = 0, left = 18, right = 21, interval = 161
fit_prep <- betaregscale(y ~ x1 + x2, data = prep)
summary(fit_prep)
#> 
#> Call:
#> betaregscale(formula = y ~ x1 + x2, data = prep)
#> 
#> Quantile residuals:
#>     Min      1Q  Median      3Q     Max 
#> -3.0625 -0.5896  0.2555  0.6723  1.5528 
#> 
#> Coefficients (mean model with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -5.15677    0.08757 -58.886  < 2e-16 ***
#> x1          -0.38266    0.06216  -6.156 7.45e-10 ***
#> x2           0.12275    0.06556   1.872   0.0612 .  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Phi coefficients (precision model with logit link):
#>       Estimate Std. Error z value Pr(>|z|)    
#> (phi)  -4.8810     0.1352  -36.11   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Log-likelihood: -896.1340 on 4 Df
#> Pseudo R-squared: 0.1097 
#> Number of iterations: 29 (BFGS) 
#> Censoring: 161 interval | 18 left | 21 right
```

## Example 1: Fixed dispersion model

### Simulating data

We simulate 200 observations from a beta regression model with fixed
dispersion, two covariates, and logit link for the mean.

``` r
set.seed(4255)
n <- 200
dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))

sim_fixed <- betaregscale_simulate(
  formula  = ~ x1 + x2,
  data     = dat,
  beta     = c(0.3, -0.6, 0.4),
  phi      = 1 / 10,
  link     = "logit",
  link_phi = "logit",
  ncuts    = 100,
  type     = "m",
  repar    = 2
)

head(sim_fixed)
#>    left   right      yt   y delta          x1          x2
#> 1 0.005 0.01500 0.01000   1     3  1.95102377 -0.55883986
#> 2 0.665 0.67500 0.67000  67     3  0.77253858 -0.07106106
#> 3 0.285 0.29500 0.29000  29     3  0.72640816  0.46916754
#> 4 0.895 0.90500 0.90000  90     3  0.04873961  0.11129810
#> 5 0.315 0.32500 0.32000  32     3 -0.54450108 -0.56115817
#> 6 0.995 0.99999 0.99999 100     2  0.36002855  0.01866439
```

The `type = "m"` argument means that each observation is centered in its
interval. For example, a score of 67 on a 0–100 scale yields
$y_{t} = 0.67$ with interval $\lbrack 0.665,0.675\rbrack$.

### Fitting the model

``` r
fit_fixed <- betaregscale(
  y ~ x1 + x2,
  data     = sim_fixed,
  link     = "logit",
  link_phi = "logit",
  repar    = 2
)
summary(fit_fixed)
#> 
#> Call:
#> betaregscale(formula = y ~ x1 + x2, data = sim_fixed, link = "logit", 
#>     link_phi = "logit", repar = 2)
#> 
#> Quantile residuals:
#>     Min      1Q  Median      3Q     Max 
#> -3.5250 -0.6758  0.0912  0.6465  2.7440 
#> 
#> Coefficients (mean model with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.33923    0.09872   3.436 0.000590 ***
#> x1          -0.64453    0.10940  -5.892 3.82e-09 ***
#> x2           0.37387    0.10820   3.455 0.000549 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Phi coefficients (precision model with logit link):
#>       Estimate Std. Error z value Pr(>|z|)  
#> (phi)  0.20493    0.09322   2.198   0.0279 *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Log-likelihood: -786.6074 on 4 Df
#> Pseudo R-squared: 0.2154 
#> Number of iterations: 40 (BFGS) 
#> Censoring: 154 interval | 13 left | 33 right
```

The summary output follows the `betareg` package style, showing separate
coefficient tables for the mean and precision submodels, with Wald
z-tests and $p$-values based on the standard normal distribution.

### Goodness of fit

``` r
gof(fit_fixed)
#>      logLik      AIC      BIC pseudo_r2
#> 1 -786.6074 1581.215 1594.408 0.2153843
```

### Comparing link functions

The package supports several link functions for the mean submodel. We
can compare them using information criteria:

``` r
links <- c("logit", "probit", "cauchit", "cloglog")
fits <- lapply(setNames(links, links), function(lnk) {
  betaregscale(y ~ x1 + x2, data = sim_fixed, link = lnk, repar = 2)
})

# Estimates
est_table <- do.call(rbind, lapply(names(fits), function(lnk) {
  e <- est(fits[[lnk]])
  e$link <- lnk
  e
}))
est_table
#>       variable   estimate         se   z_value      p_value     ci_lower
#> 1  (Intercept)  0.3392311 0.09872241  3.436212 5.899089e-04  0.145738775
#> 2           x1 -0.6445300 0.10939654 -5.891686 3.822756e-09 -0.858943312
#> 3           x2  0.3738654 0.10819762  3.455394 5.494894e-04  0.161801954
#> 4        (phi)  0.2049261 0.09322014  2.198303 2.792749e-02  0.022218033
#> 5  (Intercept)  0.2041063 0.06006633  3.398015 6.787662e-04  0.086378463
#> 6           x1 -0.3821542 0.06311114 -6.055258 1.401925e-09 -0.505849793
#> 7           x2  0.2235686 0.06498831  3.440136 5.814209e-04  0.096193905
#> 8        (phi)  0.2115155 0.09287685  2.277376 2.276375e-02  0.029480271
#> 9  (Intercept)  0.3339365 0.09607027  3.475961 5.090262e-04  0.145642243
#> 10          x1 -0.6714835 0.13361627 -5.025462 5.022215e-07 -0.933366566
#> 11          x2  0.3654699 0.10750561  3.399543 6.749860e-04  0.154762806
#> 12       (phi)  0.1830339 0.09397476  1.947693 5.145175e-02 -0.001153193
#> 13 (Intercept) -0.1693231 0.06616497 -2.559105 1.049421e-02 -0.299004071
#> 14          x1 -0.3757247 0.06296732 -5.966979 2.416865e-09 -0.499138330
#> 15          x2  0.2136847 0.06629000  3.223482 1.266421e-03  0.083758640
#> 16       (phi)  0.2323318 0.09199950  2.525359 1.155802e-02  0.052016049
#>       ci_upper    link
#> 1   0.53272352   logit
#> 2  -0.43011676   logit
#> 3   0.58592882   logit
#> 4   0.38763425   logit
#> 5   0.32183414  probit
#> 6  -0.25845868  probit
#> 7   0.35094339  probit
#> 8   0.39355082  probit
#> 9   0.52223077 cauchit
#> 10 -0.40960043 cauchit
#> 11  0.57617705 cauchit
#> 12  0.36722108 cauchit
#> 13 -0.03964214 cloglog
#> 14 -0.25231098 cloglog
#> 15  0.34361067 cloglog
#> 16  0.41264747 cloglog

# Goodness of fit
gof_table <- do.call(rbind, lapply(fits, gof))
gof_table
#>            logLik      AIC      BIC pseudo_r2
#> logit   -786.6074 1581.215 1594.408 0.2153843
#> probit  -786.9783 1581.957 1595.150 0.2137397
#> cauchit -785.1898 1578.380 1591.573 0.1729558
#> cloglog -788.2910 1584.582 1597.775 0.1589974
```

### Residual diagnostics

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method
provides six diagnostic panels. By default, the first four are shown:

``` r
plot(fit_fixed)
```

![](betaregscale_files/figure-html/plot-fixed-1.png)

For ggplot2 output (requires the **ggplot2** package):

``` r
plot(fit_fixed, gg = TRUE)
```

![](betaregscale_files/figure-html/plot-fixed-gg-1.png)

### Predictions

``` r
# Fitted means
head(predict(fit_fixed, type = "response"))
#> [1] 0.2446795 0.4538169 0.5116091 0.5864787 0.6178378 0.5285090

# Conditional variance
head(predict(fit_fixed, type = "variance"))
#> [1] 0.1018409 0.1365879 0.1376890 0.1336422 0.1301115 0.1373154

# Quantile predictions
head(predict(fit_fixed, type = "quantile", at = c(0.25, 0.5, 0.75)))
#>           q_0.25      q_0.5    q_0.75
#> [1,] 0.002041529 0.06473391 0.4288653
#> [2,] 0.073601085 0.40611615 0.8351903
#> [3,] 0.125179599 0.52371492 0.8976529
#> [4,] 0.220087703 0.67348688 0.9525970
#> [5,] 0.270331988 0.73257944 0.9680784
#> [6,] 0.143685161 0.55814195 0.9124645
```

### Confidence intervals

Wald confidence intervals based on the asymptotic normal approximation:

``` r
confint(fit_fixed)
#>                   2.5 %     97.5 %
#> (Intercept)  0.14573878  0.5327235
#> x1          -0.85894331 -0.4301168
#> x2           0.16180195  0.5859288
#> (phi)        0.02221803  0.3876343
confint(fit_fixed, model = "mean")
#>                  2.5 %     97.5 %
#> (Intercept)  0.1457388  0.5327235
#> x1          -0.8589433 -0.4301168
#> x2           0.1618020  0.5859288
```

### Censoring structure

The
[`censoring_summary()`](https://evandeilton.github.io/betaregscale/reference/censoring_summary.md)
function provides a visual and tabular overview of the censoring types
in the fitted model:

``` r
censoring_summary(fit_fixed)
```

![](betaregscale_files/figure-html/censoring-summary-1.png)

## Example 2: Variable dispersion model

In many applications, the dispersion parameter $\phi$ may depend on
covariates. The package supports variable-dispersion models using the
`Formula` package notation: `y ~ x1 + x2 | z1 + z2`, where the terms
after `|` define the linear predictor for $\phi$.

### Simulating data

``` r
set.seed(2222)
n <- 200
dat_z <- data.frame(
  x1 = rnorm(n),
  x2 = rnorm(n),
  x3 = rbinom(n, size = 1, prob = 0.5),
  z1 = rnorm(n),
  z2 = rnorm(n)
)

sim_var <- betaregscale_simulate_z(
  formula_x = ~ x1 + x2 + x3,
  formula_z = ~ z1 + z2,
  data = dat_z,
  beta = c(0.2, -0.6, 0.2, 0.2),
  zeta = c(0.2, -0.8, 0.6),
  link = "logit",
  link_phi = "logit",
  ncuts = 100,
  type = "m",
  repar = 2
)

head(sim_var)
#>    left   right      yt   y delta         x1         x2 x3          z1
#> 1 0.195 0.20500 0.20000  20     3 -0.3380621 -0.4210442  1  0.44649157
#> 2 0.015 0.02500 0.02000   2     3  0.9391643  0.6791775  0 -0.25705190
#> 3 0.325 0.33500 0.33000  33     3  1.7377190 -0.9225684  0  2.42959414
#> 4 0.725 0.73500 0.73000  73     3  0.6963261 -0.1521523  0 -0.03793774
#> 5 0.425 0.43500 0.43000  43     3  0.4622959 -0.6422296  0  0.31676341
#> 6 0.995 0.99999 0.99999 100     2 -0.3150868  0.7921496  0 -1.13668590
#>           z2
#> 1  0.6886654
#> 2  0.7827447
#> 3  1.1176679
#> 4  0.3732239
#> 5 -1.1603812
#> 6 -0.8365901
```

### Fitting the model

``` r
fit_var <- betaregscale(
  y ~ x1 + x2 | z1,
  data     = sim_var,
  link     = "logit",
  link_phi = "logit",
  repar    = 2
)
summary(fit_var)
#> 
#> Call:
#> betaregscale(formula = y ~ x1 + x2 | z1, data = sim_var, link = "logit", 
#>     link_phi = "logit", repar = 2)
#> 
#> Quantile residuals:
#>     Min      1Q  Median      3Q     Max 
#> -3.3476 -0.6251 -0.0383  0.6391  2.4092 
#> 
#> Coefficients (mean model with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.27958    0.08675   3.223  0.00127 ** 
#> x1          -0.44820    0.09340  -4.799  1.6e-06 ***
#> x2           0.19874    0.08391   2.368  0.01787 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Phi coefficients (precision model with logit link):
#>                   Estimate Std. Error z value Pr(>|z|)    
#> (phi)_(Intercept)  0.08143    0.09412   0.865    0.387    
#> (phi)_z1          -0.91702    0.10308  -8.896   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Log-likelihood: -810.2369 on 5 Df
#> Pseudo R-squared: 0.0690 
#> Number of iterations: 38 (BFGS) 
#> Censoring: 161 interval | 15 left | 24 right
```

Notice the `(phi)_` prefix in the precision coefficient names, following
the `betareg` convention.

### Accessing coefficients by submodel

``` r
# Full parameter vector
coef(fit_var)
#>       (Intercept)                x1                x2 (phi)_(Intercept) 
#>        0.27958256       -0.44820458        0.19873548        0.08142918 
#>          (phi)_z1 
#>       -0.91701510

# Mean submodel only
coef(fit_var, model = "mean")
#> (Intercept)          x1          x2 
#>   0.2795826  -0.4482046   0.1987355

# Precision submodel only
coef(fit_var, model = "precision")
#> (phi)_(Intercept)          (phi)_z1 
#>        0.08142918       -0.91701510

# Variance-covariance matrix for the mean submodel
vcov(fit_var, model = "mean")
#>               (Intercept)            x1           x2
#> (Intercept)  7.526312e-03 -0.0010535138 9.147362e-06
#> x1          -1.053514e-03  0.0087243197 1.616653e-04
#> x2           9.147362e-06  0.0001616653 7.041609e-03
```

### Comparing link functions (variable dispersion)

``` r
links <- c("logit", "probit", "cauchit", "cloglog")
fits_var <- lapply(setNames(links, links), function(lnk) {
  betaregscale(y ~ x1 + x2 | z1, data = sim_var, link = lnk, repar = 2)
})

# Estimates
est_var <- do.call(rbind, lapply(names(fits_var), function(lnk) {
  e <- est(fits_var[[lnk]])
  e$link <- lnk
  e
}))
est_var
#>             variable    estimate         se    z_value      p_value    ci_lower
#> 1        (Intercept)  0.27958256 0.08675432  3.2226932 1.269915e-03  0.10954722
#> 2                 x1 -0.44820458 0.09340407 -4.7985553 1.598142e-06 -0.63127318
#> 3                 x2  0.19873548 0.08391430  2.3683150 1.786932e-02  0.03426649
#> 4  (phi)_(Intercept)  0.08142918 0.09412224  0.8651427 3.869604e-01 -0.10304703
#> 5           (phi)_z1 -0.91701510 0.10308437 -8.8957729 5.801413e-19 -1.11905675
#> 6        (Intercept)  0.17155835 0.05343984  3.2103083 1.325927e-03  0.06681819
#> 7                 x1 -0.27344889 0.05634629 -4.8530059 1.216041e-06 -0.38388559
#> 8                 x2  0.12201196 0.05157291  2.3658148 1.799044e-02  0.02093091
#> 9  (phi)_(Intercept)  0.08344785 0.09404841  0.8872862 3.749248e-01 -0.10088364
#> 10          (phi)_z1 -0.91678540 0.10310939 -8.8913864 6.035111e-19 -1.11887609
#> 11       (Intercept)  0.24724531 0.07689481  3.2153707 1.302762e-03  0.09653426
#> 12                x1 -0.40391079 0.08959490 -4.5081894 6.538319e-06 -0.57951357
#> 13                x2  0.17683321 0.07462238  2.3697075 1.780216e-02  0.03057604
#> 14 (phi)_(Intercept)  0.07330179 0.09429400  0.7773749 4.369376e-01 -0.11151105
#> 15          (phi)_z1 -0.91651393 0.10304420 -8.8943768 5.874807e-19 -1.11847685
#> 16       (Intercept) -0.19102845 0.05819476 -3.2825715 1.028649e-03 -0.30508808
#> 17                x1 -0.28317419 0.05923492 -4.7805283 1.748351e-06 -0.39927249
#> 18                x2  0.12614368 0.05504238  2.2917555 2.191975e-02  0.01826260
#> 19 (phi)_(Intercept)  0.09029228 0.09385817  0.9620076 3.360458e-01 -0.09366636
#> 20          (phi)_z1 -0.91498567 0.10292470 -8.8898551 6.118865e-19 -1.11671437
#>       ci_upper    link
#> 1   0.44961791   logit
#> 2  -0.26513597   logit
#> 3   0.36320448   logit
#> 4   0.26590539   logit
#> 5  -0.71497346   logit
#> 6   0.27629850  probit
#> 7  -0.16301219  probit
#> 8   0.22309301  probit
#> 9   0.26777934  probit
#> 10 -0.71469472  probit
#> 11  0.39795636 cauchit
#> 12 -0.22830801 cauchit
#> 13  0.32309038 cauchit
#> 14  0.25811462 cauchit
#> 15 -0.71455101 cauchit
#> 16 -0.07696882 cloglog
#> 17 -0.16707589 cloglog
#> 18  0.23402477 cloglog
#> 19  0.27425092 cloglog
#> 20 -0.71325696 cloglog

# Goodness of fit
gof_var <- do.call(rbind, lapply(fits_var, gof))
gof_var
#>            logLik      AIC      BIC  pseudo_r2
#> logit   -810.2369 1630.474 1646.965 0.06902340
#> probit  -810.3704 1630.741 1647.232 0.08028804
#> cauchit -809.6065 1629.213 1645.705 0.03799137
#> cloglog -810.9149 1631.830 1648.321 0.03641597
```

### Diagnostics for variable dispersion

``` r
plot(fit_var)
```

![](betaregscale_files/figure-html/plot-variable-1.png)

## S3 methods reference

The following standard S3 methods are available for objects of class
`"betaregscale"`:

| Method                                                                                   | Description                                                |
|:-----------------------------------------------------------------------------------------|:-----------------------------------------------------------|
| [`print()`](https://rdrr.io/r/base/print.html)                                           | Compact display of call and coefficients                   |
| [`summary()`](https://rdrr.io/r/base/summary.html)                                       | Detailed output with Wald tests and goodness-of-fit        |
| `coef(model=)`                                                                           | Extract coefficients (full, mean, or precision)            |
| `vcov(model=)`                                                                           | Variance-covariance matrix (full, mean, or precision)      |
| `confint(model=)`                                                                        | Wald confidence intervals                                  |
| [`logLik()`](https://rdrr.io/r/stats/logLik.html)                                        | Log-likelihood value                                       |
| [`AIC()`](https://rdrr.io/r/stats/AIC.html), [`BIC()`](https://rdrr.io/r/stats/AIC.html) | Information criteria                                       |
| [`nobs()`](https://rdrr.io/r/stats/nobs.html)                                            | Number of observations                                     |
| [`formula()`](https://rdrr.io/r/stats/formula.html)                                      | Model formula                                              |
| `model.matrix(model=)`                                                                   | Design matrix (mean or precision)                          |
| [`fitted()`](https://rdrr.io/r/stats/fitted.values.html)                                 | Fitted mean values                                         |
| `residuals(type=)`                                                                       | Residuals: response, pearson, rqr, weighted, sweighted     |
| `predict(type=)`                                                                         | Predictions: response, link, precision, variance, quantile |
| `plot(gg=)`                                                                              | Diagnostic plots (base R or ggplot2)                       |

## Reparameterizations

The package supports three reparameterizations of the beta distribution,
controlled by the `repar` argument:

**Direct (`repar = 0`):** Shape parameters $a = \mu$ and $b = \phi$ are
used directly. This is rarely used in practice.

**Precision (`repar = 1`, Ferrari & Cribari-Neto, 2004):** The mean
$\mu \in (0,1)$ and precision $\phi > 0$ yield $a = \mu\phi$ and
$b = (1 - \mu)\phi$. Higher $\phi$ means less variability.

**Mean–variance (`repar = 2`):** The mean $\mu \in (0,1)$ and dispersion
$\phi \in (0,1)$ yield $a = \mu(1 - \phi)/\phi$ and
$b = (1 - \mu)(1 - \phi)/\phi$. Here $\phi$ acts as a coefficient of
variation: smaller $\phi$ means less variability.

``` r
# Precision parameterization: mu = 0.5, phi = 10 (high precision)
beta_reparam(mu = 0.5, phi = 10, repar = 1)
#>   shape1 shape2
#> 1      5      5

# Mean-variance parameterization: mu = 0.5, phi = 0.1 (low dispersion)
beta_reparam(mu = 0.5, phi = 0.1, repar = 2)
#>   shape1 shape2
#> 1    4.5    4.5
```

## References

- Lopes, J. E. (2024). *Beta Regression for Interval-Censored
  Scale-Derived Outcomes*. MSc Dissertation, PPGMNE/UFPR.

- Ferrari, S. and Cribari-Neto, F. (2004). Beta regression for modelling
  rates and proportions. *Journal of Applied Statistics*, **31**(7),
  799–815.
