# betaregscale <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen)](https://github.com/evandeilton/betaregscale)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**Beta Regression for Interval-Censored Scale-Derived Outcomes**

`betaregscale` provides maximum-likelihood estimation of beta regression
models for responses derived from bounded rating scales. Observations are
treated as interval-censored on (0, 1) and the likelihood is built from
the difference of the beta CDF at the interval endpoints.

## Key features

- **Fixed and variable dispersion** — model a scalar φ or let it depend
  on covariates via a second linear predictor (`y ~ x1 + x2 | z1`).
- **High-performance C++ backend** — the log-likelihood core and its
  gradient are compiled via Rcpp/RcppArmadillo with full numerical
  safeguards.
- **Flexible link functions** — logit, probit, cauchit, cloglog for the
  mean; logit, log, identity, sqrt, inverse, and more for φ.
- **Three reparameterizations** — direct (0), Ferrari–Cribari-Neto
  precision (1), or mean–variance (2).
- **Full S3 interface** — `print`, `summary`, `coef`, `vcov`, `logLik`,
  `AIC`, `BIC`, `fitted`, `residuals`, `predict`.
- **Simulation toolkit** — `betaregscale_simulate()` and
  `betaregscale_simulate_z()` for Monte Carlo studies.
- **bbmle integration** — profile-likelihood CIs and LR tests via
  `betaregscale_bbmle()`.

## Installation

```r
# Development version from GitHub:
# install.packages("remotes")
remotes::install_github("evandeilton/betaregscale")
```

## Quick start

```r
library(betaregscale)

# Simulate interval-censored data
set.seed(42)
n <- 500
dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
sim <- betaregscale_simulate(
  formula = ~ x1 + x2, data = dat,
  beta = c(0.2, -0.5, 0.3), phi = 1/5,
  link = "logit", link_phi = "logit"
)

# Fit the model
fit <- betaregscale(y ~ x1 + x2, data = sim)
print(fit)
summary(fit)

# Predictions
predict(fit, type = "response")
```

### Variable dispersion

```r
dat$z1 <- runif(n)
sim_z <- betaregscale_simulate_z(
  formula_x = ~ x1 + x2, formula_z = ~ z1,
  data = dat,
  beta = c(0.2, -0.5, 0.3),
  zeta = c(1, 1.2)
)
fit_z <- betaregscale(y ~ x1 + x2 | z1, data = sim_z)
print(fit_z)
```

## Model details

The interval log-likelihood for observation *i* is:

$$
\ell_i = \log\!\bigl[F_\text{Beta}(r_i;\, a_i, b_i) - F_\text{Beta}(l_i;\, a_i, b_i)\bigr]
$$

where $[l_i, r_i]$ are the interval endpoints and $(a_i, b_i)$ are the
beta shape parameters obtained from $\mu_i = g^{-1}(x_i'\beta)$ and
$\phi_i = h^{-1}(z_i'\gamma)$ via the chosen reparameterization.

## References

- Lopes, J. E. (2024). *Beta Regression for Interval-Censored
  Scale-Derived Outcomes*. MSc Dissertation, PPGMNE/UFPR.
- Ferrari, S. and Cribari-Neto, F. (2004). Beta regression for
  modelling rates and proportions. *Journal of Applied Statistics*,
  31(7), 799–815.

## License

MIT © Jose Eduardo Lopes
