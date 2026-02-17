# betaregscale <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/evandeilton/betaregscale/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/evandeilton/betaregscale/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![CRAN status](https://www.r-pkg.org/badges/version/betaregscale)](https://CRAN.R-project.org/package=betaregscale) 
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/betaregscale)](https://cran.r-project.org/package=betaregscale)
<!-- badges: end -->

**Beta Regression for Interval-Censored Scale-Derived Outcomes**

The `betaregscale` package provides maximum-likelihood estimation of beta
regression models for responses derived from bounded rating scales (e.g.,
NRS-11, NRS-21, NRS-101, Likert scales). Observations are treated as
interval-censored on (0, 1) after a scale-to-unit transformation, and
the likelihood is built from the beta CDF evaluated at the interval
endpoints.

The package is designed for situations where the recorded score carries
measurement uncertainty inherent to the instrument. For example, a pain
score of 6 on a (0-10) NRS scale is interpreted as lying in the interval
[5.5, 6.5] after rescaling to (0, 1). The complete likelihood (Lopes,
2024, Eq. 2.24) supports mixed censoring types: **uncensored,
left-censored, right-censored, and interval-censored** within the same
dataset.

Mathematically, for each observation $i$, the likelihood contribution is

$$
L_i(\theta)=
\begin{cases}
f(y_i; a_i, b_i), & \delta_i=0,\\
F(u_i; a_i, b_i), & \delta_i=1,\\
1-F(l_i; a_i, b_i), & \delta_i=2,\\
F(u_i; a_i, b_i)-F(l_i; a_i, b_i), & \delta_i=3.
\end{cases}
$$

where $f(\cdot)$ and $F(\cdot)$ denote the beta density and CDF.

## Key features

- **Mixed censoring support**: the complete likelihood handles four
  censoring types simultaneously: exact observations ($\delta=0$),
  left-censored ($\delta=1$), right-censored ($\delta=2$), and
  interval-censored ($\delta=3$).
- **Fixed and variable dispersion**: model a scalar $\phi$ or let it
  depend on covariates via a second linear predictor (`y ~ x1 + x2 | z1`).
- **Mixed-effects support**: `brsmm()` fits Gaussian mixed beta interval
  models with random effects in the mean predictor (e.g.
  `random = ~ 1 | group` and `random = ~ 1 + x | group`) via
  Laplace-approximated marginal likelihood.
- **High-performance C++ backend**: the log-likelihood and numerical
  gradient are compiled via Rcpp/RcppArmadillo for fast, numerically
  stable estimation.
- **Flexible link functions**: logit, probit, cauchit, cloglog for
  the mean; logit, log, identity, sqrt, inverse, and others for the
  dispersion.
- **Three reparameterizations**: direct (0), Ferrari--Cribari-Neto
  precision (1), and mean--variance (2).
- **betareg-style S3 interface**: `print`, `summary`, `coef`, `vcov`,
  `logLik`, `AIC`, `BIC`, `nobs`, `formula`, `model.matrix`, `fitted`,
  `residuals`, `predict`, `confint`, and `plot` methods. The `coef()`
  and `vcov()` methods accept `model = c("full", "mean", "precision")`.
- **Diagnostic plots**: six residual diagnostic panels with both base
  R and ggplot2 backends, including a half-normal envelope plot.
- **Censoring summary**: `brs_cens()` provides visual and
  tabular summaries of the censoring structure.
- **Simulation toolkit**: `brs_sim()` supports both fixed- and
  variable-dispersion Monte Carlo studies via one- or two-part formulas.
- **Analyst toolkit**: `brs_bootstrap()`, `brs_marginaleffects()`,
  `brs_predict_scoreprob()`, `brs_cv()`, and `brs_table()` for
  uncertainty quantification, interpretation, and model validation.

## Installation

```r
# Development version from GitHub:
# install.packages("remotes")
remotes::install_github("evandeilton/betaregscale")
```

## Quick start

### Fixed dispersion model

```r
library(betaregscale)

# Simulate interval-censored data from a fixed-dispersion model
set.seed(42)
n <- 200
dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
sim <- brs_sim(
  formula = ~ x1 + x2, data = dat,
  beta = c(0.3, -0.6, 0.4), phi = 1/10,
  link = "logit", link_phi = "logit",
  ncuts = 100, repar = 2
)
head(sim)

# Fit the model
fit <- brs(y ~ x1 + x2, data = sim, repar = 2)
summary(fit)
```

### Variable dispersion model

```r
# Simulate data with covariate-dependent dispersion
set.seed(2222)
n <- 200
dat <- data.frame(
  x1 = rnorm(n), x2 = rnorm(n),
  z1 = rnorm(n), z2 = rnorm(n)
)
sim_z <- brs_sim(
  formula = ~ x1 + x2 | z1,
  data = dat,
  beta = c(0.2, -0.6, 0.2),
  zeta = c(0.2, -0.8),
  link = "logit", link_phi = "logit",
  ncuts = 100, repar = 2
)

# Fit variable-dispersion model (pipe notation)
fit_z <- brs(y ~ x1 + x2 | z1, data = sim_z, repar = 2)
summary(fit_z)
```

### Comparing link functions

```r
links <- c("logit", "probit", "cauchit", "cloglog")
fits <- lapply(setNames(links, links), function(lnk) {
  brs(y ~ x1 + x2, data = sim, link = lnk, repar = 2)
})

# Goodness-of-fit comparison
do.call(rbind, lapply(fits, brs_gof))

# Analyst-friendly table
knitr::kable(do.call(rbind, lapply(fits, brs_gof)), digits = 4)
```

### Mixed model (random intercept)

```r
set.seed(123)
g <- 12
ni <- 8
id <- factor(rep(seq_len(g), each = ni))
n <- length(id)
x1 <- rnorm(n)
b <- rnorm(g, sd = 0.5)
mu <- plogis(0.2 + 0.6 * x1 + b[as.integer(id)])
phi <- plogis(-0.3 + 0.2 * x1)
shp <- brs_repar(mu = mu, phi = phi, repar = 2)
y <- round(rbeta(n, shp$shape1, shp$shape2) * 100)
dmm <- data.frame(y = y, x1 = x1, id = id)

fit_mm <- brsmm(y ~ x1, random = ~ 1 | id, data = dmm, repar = 2)
summary(fit_mm)
```

### Mixed model (random intercept + random slope)

```r
fit_mm_rs <- brsmm(
  y ~ x1,
  random = ~ 1 + x1 | id,
  data = dmm,
  repar = 2,
  int_method = "laplace"
)

summary(fit_mm_rs)
knitr::kable(fit_mm_rs$random$D, digits = 4)
knitr::kable(head(ranef.brsmm(fit_mm_rs)), digits = 4)
```

### Evolutionary model selection (`anova`)

```r
fit_brs <- brs(y ~ x1, data = dmm, repar = 2)
fit_mm_ri <- brsmm(y ~ x1, random = ~ 1 | id, data = dmm, repar = 2)
fit_mm_rs <- brsmm(y ~ x1, random = ~ 1 + x1 | id, data = dmm, repar = 2)

anova(fit_brs, fit_mm_ri, fit_mm_rs, test = "Chisq")
```

### Random-effects numerical study

```r
re_study <- brsmm_re_study(fit_mm_rs)
print(re_study)
knitr::kable(re_study$summary, digits = 4)
```

### S3 methods

```r
# Coefficients by submodel
coef(fit)                          # full parameter vector
coef(fit, model = "mean")         # mean submodel only
coef(fit, model = "precision")    # precision submodel only

# Variance-covariance matrix
vcov(fit, model = "mean")

# Wald confidence intervals
confint(fit)

# Predictions
predict(fit, type = "response")    # fitted means
predict(fit, type = "precision")   # fitted dispersion
predict(fit, type = "variance")    # conditional variance
predict(fit, type = "quantile", at = c(0.25, 0.5, 0.75))

# Residuals
residuals(fit, type = "pearson")
residuals(fit, type = "rqr")       # randomized quantile residuals

# Diagnostic plots (base R)
plot(fit)

# Diagnostic plots (ggplot2, if installed)
plot(fit, gg = TRUE)

# Censoring structure summary
brs_cens(fit)
```

## Model details

### Complete likelihood

The complete log-likelihood for mixed censoring (Lopes, 2024, Eq. 2.24)
combines the four observation types:

$$
\ell(\theta)=
\sum_{i:\delta_i=0}\log f(y_i)+
\sum_{i:\delta_i=1}\log F(u_i)+
\sum_{i:\delta_i=2}\log\!\left[1-F(l_i)\right]+
\sum_{i:\delta_i=3}\log\!\left[F(u_i)-F(l_i)\right].
$$

where $f(\cdot)$ and $F(\cdot)$ are the beta density and CDF,
$[l_i, u_i]$ are the interval endpoints, and $\delta_i$ indicates the
censoring type.

### Reparameterizations

| Code | Name | Shape parameters |
|------|------|-----------------|
| 0 | Direct | $a = \mu,\; b = \phi$ |
| 1 | Precision (Ferrari & Cribari-Neto) | $a = \mu\phi,\; b = (1-\mu)\phi$ |
| 2 | Mean--variance | $a = \mu(1-\phi)/\phi,\; b = (1-\mu)(1-\phi)/\phi$ |

### Analyst-oriented outputs (clean tables)

```r
# Parametric bootstrap CI
boot_ci <- brs_bootstrap(fit, R = 100, level = 0.95)
knitr::kable(head(boot_ci), digits = 4)

# Average marginal effects
ame <- brs_marginaleffects(fit, model = "mean", type = "response")
knitr::kable(ame, digits = 4)

# Score-scale probabilities
ps <- brs_predict_scoreprob(fit, scores = 0:10)
knitr::kable(ps[1:6, 1:6], digits = 4)

# Cross-validation summary
cv <- brs_cv(y ~ x1 + x2, data = sim, k = 5, repeats = 1, seed = 123)
knitr::kable(cv, digits = 4)
```

### Interval construction

Scale observations are mapped to (0, 1) with midpoint uncertainty intervals:

$$
y_t = y/K, \quad \text{interval } [y_t - h/K,\; y_t + h/K]
$$

where $K$ is the number of scale categories (`ncuts`) and $h$ is the
half-width (`lim`, default `0.5`).

### Mixed-model likelihood (`brsmm`)

For group $j$ with random-effects vector $\mathbf{b}_j \in \mathbb{R}^{q_b}$:

$$
\eta_{\mu,ij} = x_{ij}^\top\beta + w_{ij}^\top\mathbf{b}_j,
\qquad
\eta_{\phi,ij} = z_{ij}^\top\gamma,
$$

with

$$
\mathbf{b}_j \sim \mathcal{N}(\mathbf{0}, D).
$$

The marginal group likelihood is

$$
L_j(\theta)=\int_{\mathbb{R}^{q_b}}
\left\{\prod_{i=1}^{n_j} L_{ij}(\mathbf{b}_j;\theta)\right\}
\varphi_{q_b}(\mathbf{b}_j;\mathbf{0},D)\,d\mathbf{b}_j,
$$

and the log-likelihood is $\ell(\theta)=\sum_{j=1}^G \log L_j(\theta)$.
`brsmm()` uses a multivariate Laplace approximation:

$$
\log L_j(\theta)\approx
Q_j(\hat{\mathbf{b}}_j)+\frac{q_b}{2}\log(2\pi)-\frac{1}{2}\log|H_j|,
$$

where
$Q_j(\mathbf{b})=\sum_i \log L_{ij}(\mathbf{b};\theta)+\log\varphi_{q_b}(\mathbf{b};\mathbf{0},D)$
and $H_j=-\nabla^2Q_j(\hat{\mathbf{b}}_j)$.

## References

- Lopes, J. E. (2024). *Beta Regression for Interval-Censored
  Scale-Derived Outcomes*. MSc Dissertation, PPGMNE/UFPR.
- Ferrari, S. L. P. and Cribari-Neto, F. (2004). Beta regression for
  modelling rates and proportions. *Journal of Applied Statistics*,
  31(7), 799--815. DOI: 10.1080/0266476042000214501.
  DOI validated via: <https://doi.org/10.1080/0266476042000214501>.

- Hawker, G. A., Mian, S., Kendzerska, T., and French, M. (2011).
  Measures of adult pain: VAS, NRS, MPQ, SF-MPQ, CPGS, SF-36 BPS,
  and ICOAP. *Arthritis Care and Research*, 63(S11), S240--S252.
  DOI: 10.1002/acr.20543. Validated online via:
  <https://doi.org/10.1002/acr.20543>.

- Hjermstad, M. J., Fayers, P. M., Haugen, D. F., et al. (2011).
  Studies comparing numerical rating scales, verbal rating scales, and
  visual analogue scales for assessment of pain intensity in adults.
  *Journal of Pain and Symptom Management*, 41(6), 1073--1093.
  DOI: 10.1016/j.jpainsymman.2010.08.016. Validated online via:
  <https://doi.org/10.1016/j.jpainsymman.2010.08.016>.

## License

MIT &copy; José Evandeilton Lopes
