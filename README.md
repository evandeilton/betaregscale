# betaregscale <img src="man/figures/logo.png" align="right" height="139" />

[![R-CMD-check](https://github.com/evandeilton/betaregscale/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/evandeilton/betaregscale/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/betaregscale)](https://CRAN.R-project.org/package=betaregscale)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

**betaregscale** fits regression models for scores on a bounded rating scale
$0, 1, \ldots, K$ (pain NRS-11 or NRS-101, Likert items, quality ratings). Each
score is the coarsened reading of a latent beta variable on $(0, 1)$, written in
mean–dispersion form: mean $\mu$, dispersion $\sigma = 1/(1 + \phi) \in (0, 1)$
with $\phi$ the beta precision, and $\mathrm{Var}(Y) = \mu(1 - \mu)\sigma$. Both
can depend on covariates through link functions. The latent scale is cut at
cutpoints anchored on the score labels, so each score is an interval-censored
observation and the likelihood is a probability mass function over the scores
$0, \ldots, K$. The scale borders need no ad hoc squeezing: score 0 is
left-censored and score $K$ right-censored. Fixed-effects models are fitted with
`brs()` and mixed models with Gaussian random effects with `brsmm()`; the
log-likelihood and its derivatives are compiled C++.

## Installation

```r
# Released version from CRAN
install.packages("betaregscale")

# Development version from GitHub
# install.packages("remotes")
remotes::install_github("evandeilton/betaregscale")
```

## Example

Simulated NRS-11 scores (0 to 10) for 200 patients, with a continuous
covariate and a treatment indicator:

```r
library(betaregscale)

set.seed(1)
d <- data.frame(x1 = rnorm(200), trt = rbinom(200, 1, 0.5))
nrs <- brs_sim(~ x1 + trt, data = d, beta = c(-0.4, 0.5, -0.6),
               phi = qlogis(0.15), ncuts = 10)

# Mean with logit link, constant dispersion with logit link
fit <- brs(y ~ x1 + trt, data = nrs)
summary(fit)

# The fitted model is a pmf over the scores 0..10
round(head(brs_predict_scoreprob(fit), 3), 3)

# Expected score by arm at x1 = 0, and the difference between arms,
# with delta-method standard errors
nd <- data.frame(x1 = 0, trt = c(0, 1))
brs_marginal(fit, nd, type = "expected_score")
brs_marginal(fit, nd, type = "expected_score", contrast = c(2, 1))

# Diagnostic plots of randomized quantile residuals
plot(fit)
```

## Features

- **Interval maps** `interval = "mid"`, `"right"` or `"left"`, and
  **estimable boundary cutpoints** (`cutpoints = "boundary"`) for scores that
  pile up at the ends of the scale, with a likelihood-ratio test against the
  fixed map.
- **Variable dispersion** with a two-part formula `y ~ x | z`, and three
  parameterizations of the beta law (`repar`).
- **Mixed models** (`brsmm()`) with random intercepts and slopes, fitted by
  Laplace approximation, adaptive Gauss–Hermite quadrature or quasi-Monte
  Carlo, with chi-bar-square likelihood-ratio tests for variance components.
- **Score-scale predictions**: score probabilities
  (`brs_predict_scoreprob()`), and score probabilities, expected scores,
  exceedance probabilities and contrasts with standard errors
  (`brs_marginal()`), population-averaged for mixed models.
- **Diagnostics and validation**: randomized quantile and cell deviance
  residuals, parametric bootstrap, k-fold and group-level cross-validation,
  and simulation (`brs_sim()`, `brsmm_sim()`, `simulate()`).

## Documentation

- [Introduction](https://evandeilton.github.io/betaregscale/articles/brs-intro.html):
  the model, interval maps, the likelihood and a first fit.
- [Mixed models](https://evandeilton.github.io/betaregscale/articles/brs-mm.html):
  `brsmm()`, integration methods, random effects and marginal predictions.
- [Analyst tools](https://evandeilton.github.io/betaregscale/articles/brs-analyst-tools.html):
  predictions, contrasts, residuals, bootstrap, cross-validation, tables and plots.
- [Advanced workflows](https://evandeilton.github.io/betaregscale/articles/brs-advanced-workflows.html):
  variable dispersion, boundary cutpoints and model comparison.

Function reference and changelog: <https://evandeilton.github.io/betaregscale/>.
Bug reports and questions: <https://github.com/evandeilton/betaregscale/issues>.

## Citation

`citation("betaregscale")` gives the reference for the package. The model is
described in

- Lopes, J. E. (2023). *Modelos de regressão beta para dados de escala*.
  Master's dissertation, Universidade Federal do Paraná, Curitiba.
  <https://hdl.handle.net/1884/86624>

Other references:

- Bayer, F. M. (2011). *Modelagem e inferência em regressão beta*. PhD thesis,
  Universidade Federal de Pernambuco.
- Ferrari, S. L. P., and Cribari-Neto, F. (2004). Beta regression for modelling
  rates and proportions. *Journal of Applied Statistics*, 31(7), 799–815.
  <https://doi.org/10.1080/0266476042000214501>

## License

MIT © José Evandeilton Lopes and Wagner Hugo Bonat
