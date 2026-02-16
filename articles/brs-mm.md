# Mixed-Effects Beta Interval Regression with brsmm

## Overview

[`brsmm()`](https://evandeilton.github.io/betaregscale/reference/brsmm.md)
extends
[`brs()`](https://evandeilton.github.io/betaregscale/reference/brs.md)
to clustered data by adding a Gaussian random intercept in the mean
submodel while preserving the interval-censored beta likelihood for
scale-derived outcomes.

This vignette covers:

1.  full model mathematics;
2.  estimation by marginal maximum likelihood (Laplace approximation);
3.  practical use of all current `brsmm` methods;
4.  inferential and validation workflows, including parameter recovery.

``` r
library(betaregscale)
```

## Mathematical model

Assume observations $i = 1,\ldots,n_{j}$ within groups $j = 1,\ldots,G$,
with random intercept
$b_{j} \sim \mathcal{N}\left( 0,\sigma_{b}^{2} \right)$.

### Linear predictors

$$\eta_{\mu,ij} = x_{ij}^{\top}\beta + b_{j},\qquad\eta_{\phi,ij} = z_{ij}^{\top}\gamma$$

$$\mu_{ij} = g^{- 1}\left( \eta_{\mu,ij} \right),\qquad\phi_{ij} = h^{- 1}\left( \eta_{\phi,ij} \right)$$

with $g( \cdot )$ and $h( \cdot )$ chosen by `link` and `link_phi`.

### Beta parameterization

For each $\left( \mu_{ij},\phi_{ij} \right)$, `repar` maps to beta shape
parameters $\left( a_{ij},b_{ij} \right)$ via
[`brs_repar()`](https://evandeilton.github.io/betaregscale/reference/brs_repar.md).

### Conditional contribution by censoring type

Each observation contributes:

$$L_{ij}\left( b_{j};\theta \right) = \begin{cases}
{f\left( y_{ij};a_{ij},b_{ij} \right),} & {\delta_{ij} = 0} \\
{F\left( u_{ij};a_{ij},b_{ij} \right),} & {\delta_{ij} = 1} \\
{1 - F\left( l_{ij};a_{ij},b_{ij} \right),} & {\delta_{ij} = 2} \\
{F\left( u_{ij};a_{ij},b_{ij} \right) - F\left( l_{ij};a_{ij},b_{ij} \right),} & {\delta_{ij} = 3}
\end{cases}$$

where $l_{ij},u_{ij}$ are interval endpoints on $(0,1)$, $f( \cdot )$ is
beta density, and $F( \cdot )$ is beta CDF.

### Group marginal likelihood

$$L_{j}(\theta) = \int_{\mathbb{R}}\prod\limits_{i = 1}^{n_{j}}L_{ij}\left( b_{j};\theta \right)\,\varphi\left( b_{j};0,\sigma_{b}^{2} \right)\, db_{j}$$

$$\ell(\theta) = \sum\limits_{j = 1}^{G}\log L_{j}(\theta)$$

### Laplace approximation used by `brsmm()`

Define
$$Q_{j}(b) = \sum\limits_{i = 1}^{n_{j}}\log L_{ij}(b;\theta) + \log\varphi\left( b;0,\sigma_{b}^{2} \right)$$
and ${\widehat{b}}_{j} = \arg\max_{b}Q_{j}(b)$, with curvature
$H_{j} = - Q_{j}''\left( {\widehat{b}}_{j} \right)$. Then

$$\log L_{j}(\theta) \approx Q_{j}\left( {\widehat{b}}_{j} \right) + \frac{1}{2}\log(2\pi) - \frac{1}{2}\log\left( H_{j} \right).$$

[`brsmm()`](https://evandeilton.github.io/betaregscale/reference/brsmm.md)
maximizes the approximated $\ell(\theta)$ with
[`stats::optim()`](https://rdrr.io/r/stats/optim.html), and computes
group-level posterior modes ${\widehat{b}}_{j}$ and local standard
deviations from the local curvature.

## Simulating clustered scale data

The next helper simulates data from a known mixed model to illustrate
fitting, inference, and recovery checks.

``` r
sim_brsmm_data <- function(seed = 3501L, g = 24L, ni = 12L,
                           beta = c(0.20, 0.65),
                           gamma = c(-0.15),
                           sigma_b = 0.55) {
  set.seed(seed)
  id <- factor(rep(seq_len(g), each = ni))
  n <- length(id)
  x1 <- rnorm(n)

  b_true <- rnorm(g, mean = 0, sd = sigma_b)
  eta_mu <- beta[1] + beta[2] * x1 + b_true[as.integer(id)]
  eta_phi <- rep(gamma[1], n)

  mu <- plogis(eta_mu)
  phi <- plogis(eta_phi)
  shp <- brs_repar(mu = mu, phi = phi, repar = 2)
  y <- round(stats::rbeta(n, shp$shape1, shp$shape2) * 100)

  list(
    data = data.frame(y = y, x1 = x1, id = id),
    truth = list(beta = beta, gamma = gamma, sigma_b = sigma_b, b = b_true)
  )
}

sim <- sim_brsmm_data()
str(sim$data)
#> 'data.frame':    288 obs. of  3 variables:
#>  $ y : num  18 29 12 43 56 55 3 21 4 4 ...
#>  $ x1: num  -0.3677 -2.0069 -0.0469 -0.2468 0.7634 ...
#>  $ id: Factor w/ 24 levels "1","2","3","4",..: 1 1 1 1 1 1 1 1 1 1 ...
```

## Fitting `brsmm()`

``` r
fit_mm <- brsmm(
  y ~ x1,
  random = ~ 1 | id,
  data = sim$data,
  repar = 2,
  int_method = "laplace",
  method = "BFGS",
  control = list(maxit = 1000)
)

fit_mm
#> 
#> Call:
#> brsmm(formula = y ~ x1, random = ~1 | id, data = sim$data, repar = 2, 
#>     int_method = "laplace", method = "BFGS", control = list(maxit = 1000))
#> 
#> Mixed beta interval model (Laplace)
#> Observations: 288  | Groups: 24 
#> Log-likelihood: -1231.3848 
#> Random SD: 0.4303 
#> Convergence code: 0
```

## Core methods

### Coefficients and random effects

`coef(fit_mm, model = "random")` is on the optimizer scale
($\log\sigma_{b}$); use [`exp()`](https://rdrr.io/r/base/Log.html) to
report $\sigma_{b}$.

``` r
coef(fit_mm, model = "full")
#>       (Intercept)                x1 (phi)_(Intercept)           (sd)_id 
#>        0.01051336        0.62964563       -0.13626654       -0.84336505
coef(fit_mm, model = "mean")
#> (Intercept)          x1 
#>  0.01051336  0.62964563
coef(fit_mm, model = "precision")
#> (phi)_(Intercept) 
#>        -0.1362665
coef(fit_mm, model = "random")
#>   (sd)_id 
#> -0.843365
exp(coef(fit_mm, model = "random"))
#>   (sd)_id 
#> 0.4302602
head(ranef.brsmm(fit_mm))
#>           1           2           3           4           5           6 
#> -0.37076507  0.45544785 -0.17418739 -0.45940282  0.02481745  0.42957923
```

### Variance-covariance, summary and likelihood criteria

``` r
vc <- vcov(fit_mm)
dim(vc)
#> [1] 4 4

sm <- summary(fit_mm)
sm
#> 
#> Call:
#> brsmm(formula = y ~ x1, random = ~1 | id, data = sim$data, repar = 2, 
#>     int_method = "laplace", method = "BFGS", control = list(maxit = 1000))
#> 
#> Mixed beta interval model (Laplace)
#> Observations: 288  | Groups: 24 
#> logLik =-1231.3848 | AIC =2470.7696 | BIC =2485.4215
#> 
#>                   Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)        0.01051    0.05998   0.175   0.8609    
#> x1                 0.62965    0.08650   7.280 3.35e-13 ***
#> (phi)_(Intercept) -0.13627    0.07704  -1.769   0.0769 .  
#> (sd)_id           -0.84337    0.26045  -3.238   0.0012 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

logLik(fit_mm)
#> 'log Lik.' -1231.385 (df=4)
AIC(fit_mm)
#> [1] 2470.77
BIC(fit_mm)
#> [1] 2485.421
nobs(fit_mm)
#> [1] 288
```

### Fitted values, prediction and residuals

``` r
head(fitted(fit_mm, type = "mu"))
#> [1] 0.3562286 0.1646701 0.4037733 0.3738710 0.5300708 0.3169583
head(fitted(fit_mm, type = "phi"))
#> [1] 0.465986 0.465986 0.465986 0.465986 0.465986 0.465986

head(predict(fit_mm, type = "response"))
#>         1         1         1         1         1         1 
#> 0.3562286 0.1646701 0.4037733 0.3738710 0.5300708 0.3169583
head(predict(fit_mm, type = "link"))
#>          1          1          1          1          1          1 
#> -0.5917708 -1.6238828 -0.3897673 -0.5156458  0.1204284 -0.7677855
head(predict(fit_mm, type = "precision"))
#> [1] 0.465986 0.465986 0.465986 0.465986 0.465986 0.465986
head(predict(fit_mm, type = "variance"))
#> [1] 0.10686447 0.06409816 0.11218166 0.10908334 0.11607513 0.10088398

head(residuals(fit_mm, type = "response"))
#> [1] -0.17622865  0.12532992 -0.28377333  0.05612904  0.02992923  0.23304166
head(residuals(fit_mm, type = "pearson"))
#> [1] -0.53908822  0.49503053 -0.84724816  0.16994500  0.08784679  0.73370665
```

### Diagnostic plotting methods

[`plot.brsmm()`](https://evandeilton.github.io/betaregscale/reference/plot.brsmm.md)
supports base and ggplot2 backends:

``` r
plot(fit_mm, which = 1:4, type = "pearson")
```

![](brs-mm_files/figure-html/methods-plot-1.png)

``` r

if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot(fit_mm, which = 1:2, gg = TRUE)
}
```

![](brs-mm_files/figure-html/methods-plot-2.png)

[`autoplot.brsmm()`](https://evandeilton.github.io/betaregscale/reference/autoplot.brsmm.md)
provides focused ggplot diagnostics:

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  autoplot.brsmm(fit_mm, type = "calibration")
  autoplot.brsmm(fit_mm, type = "score_dist")
  autoplot.brsmm(fit_mm, type = "ranef_qq")
  autoplot.brsmm(fit_mm, type = "residuals_by_group")
}
```

![](brs-mm_files/figure-html/methods-autoplot-1.png)

### Prediction with `newdata`

If `newdata` contains unseen groups,
[`predict.brsmm()`](https://evandeilton.github.io/betaregscale/reference/predict.brsmm.md)
uses random effect equal to zero for those levels.

``` r
nd <- sim$data[1:8, c("x1", "id")]
predict(fit_mm, newdata = nd, type = "response")
#> [1] 0.3562286 0.1646701 0.4037733 0.3738710 0.5300708 0.3169583 0.4348323
#> [8] 0.3063863

nd_unseen <- nd
nd_unseen$id <- factor(rep("new_cluster", nrow(nd_unseen)))
predict(fit_mm, newdata = nd_unseen, type = "response")
#> [1] 0.4449724 0.2221609 0.4952496 0.4638431 0.6203876 0.4020284 0.5271241
#> [8] 0.3902401
```

## Statistical tests and validation workflow

### Wald tests (from `summary`)

[`summary.brsmm()`](https://evandeilton.github.io/betaregscale/reference/summary.brsmm.md)
reports Wald $z$-tests for each parameter:
$$z_{k} = {\widehat{\theta}}_{k}/{SE}\left( {\widehat{\theta}}_{k} \right).$$

``` r
summary(fit_mm)$coefficients
#>                      Estimate Std. Error    z value     Pr(>|z|)
#> (Intercept)        0.01051336 0.05997713  0.1752895 8.608522e-01
#> x1                 0.62964563 0.08649567  7.2795044 3.350653e-13
#> (phi)_(Intercept) -0.13626654 0.07704204 -1.7687296 7.693901e-02
#> (sd)_id           -0.84336505 0.26045257 -3.2380753 1.203390e-03
```

### Likelihood-ratio test for nested fixed effects

For nested models with the same random-effect structure, an LR test can
be computed as:
$$LR = 2\{\ell\left( {\widehat{\theta}}_{1} \right) - \ell\left( {\widehat{\theta}}_{0} \right)\},$$
using a $\chi^{2}$ reference with parameter-count difference.

``` r
fit_null <- brsmm(
  y ~ 1,
  random = ~ 1 | id,
  data = sim$data,
  repar = 2,
  int_method = "laplace",
  method = "BFGS",
  control = list(maxit = 1000)
)

ll0 <- logLik(fit_null)
ll1 <- logLik(fit_mm)
lr_stat <- 2 * (as.numeric(ll1) - as.numeric(ll0))
df_diff <- attr(ll1, "df") - attr(ll0, "df")
p_value <- stats::pchisq(lr_stat, df = df_diff, lower.tail = FALSE)

c(LR = lr_stat, df = df_diff, p_value = p_value)
```

### Residual diagnostics (quick checks)

``` r
r <- residuals(fit_mm, type = "pearson")
c(
  mean = mean(r),
  sd = stats::sd(r),
  q025 = as.numeric(stats::quantile(r, 0.025)),
  q975 = as.numeric(stats::quantile(r, 0.975))
)
#>        mean          sd        q025        q975 
#>  0.03803694  0.96104685 -1.56673066  1.62863050
```

## Parameter recovery experiment

A single-fit recovery table can be produced directly from the previous
fit:

``` r
est <- c(
  beta0 = coef(fit_mm, model = "mean")[1],
  beta1 = coef(fit_mm, model = "mean")[2],
  sigma_b = exp(coef(fit_mm, model = "random"))
)

true <- c(
  beta0 = sim$truth$beta[1],
  beta1 = sim$truth$beta[2],
  sigma_b = sim$truth$sigma_b
)

recovery_table <- data.frame(
  parameter = names(true),
  true = as.numeric(true),
  estimate = as.numeric(est[names(true)]),
  bias = as.numeric(est[names(true)] - true)
)
recovery_table
#>   parameter true estimate bias
#> 1     beta0 0.20       NA   NA
#> 2     beta1 0.65       NA   NA
#> 3   sigma_b 0.55       NA   NA
```

For a Monte Carlo recovery study, repeat simulation and fitting across
replicates:

``` r
mc_recovery <- function(R = 50L, seed = 7001L) {
  set.seed(seed)
  out <- vector("list", R)

  for (r in seq_len(R)) {
    sim_r <- sim_brsmm_data(seed = seed + r)
    fit_r <- brsmm(
      y ~ x1,
      random = ~ 1 | id,
      data = sim_r$data,
      repar = 2,
      int_method = "laplace",
      method = "BFGS",
      control = list(maxit = 1000)
    )

    out[[r]] <- c(
      beta0 = coef(fit_r, model = "mean")[1],
      beta1 = coef(fit_r, model = "mean")[2],
      sigma_b = exp(coef(fit_r, model = "random"))
    )
  }

  est <- do.call(rbind, out)
  truth <- c(beta0 = 0.20, beta1 = 0.65, sigma_b = 0.55)

  data.frame(
    parameter = colnames(est),
    truth = as.numeric(truth[colnames(est)]),
    mean_est = colMeans(est),
    bias = colMeans(est) - truth[colnames(est)],
    rmse = sqrt(colMeans((sweep(est, 2, truth[colnames(est)], "-"))^2))
  )
}

mc_recovery(R = 50)
```

## How this maps to automated package tests

The package test suite includes dedicated `brsmm` tests for:

1.  fitting with Laplace integration;
2.  one- and two-part formulas;
3.  S3 methods (`coef`, `vcov`, `summary`, `predict`, `residuals`,
    `ranef`);
4.  parameter recovery under known DGP settings.

Run locally:

``` r
devtools::test(filter = "brsmm")
```

## References

Ferrari, S. L. P. and Cribari-Neto, F. (2004). Beta regression for
modelling rates and proportions. *Journal of Applied Statistics*, 31(7),
799-815.

Lopes, J. E. (2024). *Beta Regression for Interval-Censored
Scale-Derived Outcomes*. MSc Dissertation, PPGMNE/UFPR.

Pinheiro, J. C. and Bates, D. M. (2000). *Mixed-Effects Models in S and
S-PLUS*. Springer.
