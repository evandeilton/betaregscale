# betaregscale: Beta Regression for Interval-Censored Scale-Derived Outcomes

Maximum-likelihood estimation of beta regression models for responses
derived from bounded rating scales. Observations are treated as
interval-censored on (0, 1) after a scale-to-unit transformation. The
complete likelihood supports mixed censoring types: uncensored (exact),
left-censored, right-censored, and interval-censored observations
(Lopes, 2024, Eq. 2.24). Both fixed- and variable-dispersion submodels
are supported, with flexible link functions for the mean and precision
components. A compiled C++ backend (via Rcpp and RcppArmadillo) provides
numerically stable, high-performance log-likelihood evaluation.

## Main functions

- [`betaregscale`](https://evandeilton.github.io/betaregscale/reference/betaregscale.md):

  Unified fitting interface for both fixed- and variable-dispersion
  models.

- [`betaregscale_fit`](https://evandeilton.github.io/betaregscale/reference/betaregscale_fit.md):

  Fit a fixed-dispersion model.

- [`betaregscale_fit_z`](https://evandeilton.github.io/betaregscale/reference/betaregscale_fit_z.md):

  Fit a variable-dispersion model.

- [`betaregscale_loglik`](https://evandeilton.github.io/betaregscale/reference/betaregscale_loglik.md):

  Compute the log-likelihood (fixed dispersion).

- [`betaregscale_loglik_z`](https://evandeilton.github.io/betaregscale/reference/betaregscale_loglik_z.md):

  Compute the log-likelihood (variable dispersion).

- [`betaregscale_simulate`](https://evandeilton.github.io/betaregscale/reference/betaregscale_simulate.md):

  Simulate interval-censored data from a fixed-dispersion beta model.

- [`betaregscale_simulate_z`](https://evandeilton.github.io/betaregscale/reference/betaregscale_simulate_z.md):

  Simulate data from a variable-dispersion beta model.

- [`censoring_summary`](https://evandeilton.github.io/betaregscale/reference/censoring_summary.md):

  Visual and tabular summary of censoring structure.

- [`bs_prepare`](https://evandeilton.github.io/betaregscale/reference/bs_prepare.md):

  Pre-process analyst data (validate, classify censoring, and rescale)
  before model fitting.

## S3 methods

Objects of class `"betaregscale"` support: `print`, `summary`, `coef`,
`vcov`, `logLik`, `AIC`, `BIC`, `nobs`, `formula`, `model.matrix`,
`fitted`, `residuals`, `predict`, `confint`, and `plot`.

The [`coef()`](https://rdrr.io/r/stats/coef.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) methods accept a
`model = c("full", "mean", "precision")` argument following the betareg
package convention.

## Censoring types

The complete likelihood (Lopes, 2024, Eq. 2.24) supports four censoring
types, classified automatically by
[`check_response`](https://evandeilton.github.io/betaregscale/reference/check_response.md):

- \\\delta = 0\\ (exact):

  Continuous observations in (0, 1).

- \\\delta = 1\\ (left-censored):

  Observations at the scale minimum (y = 0).

- \\\delta = 2\\ (right-censored):

  Observations at the scale maximum (y = ncuts).

- \\\delta = 3\\ (interval-censored):

  Standard scale observations between the boundaries.

## References

Lopes, J. E. (2024). *Beta Regression for Interval-Censored
Scale-Derived Outcomes*. MSc Dissertation, PPGMNE/UFPR.

Ferrari, S. and Cribari-Neto, F. (2004). Beta regression for modelling
rates and proportions. *Journal of Applied Statistics*, **31**(7),
799–815.

## See also

Useful links:

- <https://evandeilton.github.io/betaregscale/>

- Report bugs at <https://github.com/evandeilton/betaregscale/issues>

## Author

Jose Eduardo Lopes <evandeilton@gmail.com>
