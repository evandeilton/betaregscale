# betaregscale: Beta Regression for Interval-Censored Scale-Derived Outcomes

Maximum-likelihood estimation of beta regression models for responses
derived from bounded rating scales. Observations are treated as
interval-censored on (0, 1) after a scale-to-unit transformation. The
complete likelihood supports mixed censoring types: uncensored (exact),
left-censored, right-censored, and interval-censored observations. Both
fixed- and variable-dispersion submodels are supported, with flexible
link functions for the mean and precision components. A compiled C++
backend (via Rcpp and RcppArmadillo) provides numerically stable,
high-performance log-likelihood evaluation. Standard S3 methods
(print(), summary(), coef(), fitted(), residuals(), predict(), plot(),
confint(), vcov(), logLik(), AIC(), BIC()) are available for fitted
objects.

## Main functions

- [`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md):

  Unified fitting interface for both fixed- and variable-dispersion
  models.

- [`brs_fit_fixed`](https://evandeilton.github.io/betaregscale/reference/brs_fit_fixed.md):

  Fit a fixed-dispersion model.

- [`brs_fit_var`](https://evandeilton.github.io/betaregscale/reference/brs_fit_var.md):

  Fit a variable-dispersion model.

- [`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md):

  Fit a mixed-effects beta interval model with Gaussian random
  intercepts.

- [`brs_sim`](https://evandeilton.github.io/betaregscale/reference/brs_sim.md):

  Simulate interval-censored data from fixed or variable-dispersion beta
  models.

- [`brs_bootstrap`](https://evandeilton.github.io/betaregscale/reference/brs_bootstrap.md):

  Parametric bootstrap confidence intervals for `brs` model parameters.

- [`brs_cens`](https://evandeilton.github.io/betaregscale/reference/brs_cens.md):

  Visual and tabular summary of censoring structure.

- [`brs_prep`](https://evandeilton.github.io/betaregscale/reference/brs_prep.md):

  Pre-process analyst data (validate, classify censoring, and rescale)
  before model fitting.

## S3 methods

Objects of class `"brs"` support:
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`AIC()`](https://rdrr.io/r/stats/AIC.html),
[`BIC()`](https://rdrr.io/r/stats/AIC.html),
[`nobs()`](https://rdrr.io/r/stats/nobs.html),
[`formula()`](https://rdrr.io/r/stats/formula.html),
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`confint()`](https://rdrr.io/r/stats/confint.html), and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

The [`coef()`](https://rdrr.io/r/stats/coef.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) methods accept a
`model = c("full", "mean", "precision")` argument following the betareg
package convention.

## Censoring types

The complete likelihood supports four censoring types, classified
automatically by
[`brs_check`](https://evandeilton.github.io/betaregscale/reference/brs_check.md):

- \\\delta = 0\\ (exact):

  Continuous observations in (0, 1).

- \\\delta = 1\\ (left-censored):

  Observations at the scale minimum (y = 0).

- \\\delta = 2\\ (right-censored):

  Observations at the scale maximum (y = ncuts).

- \\\delta = 3\\ (interval-censored):

  Standard scale observations between the boundaries.

## References

Lopes, J. E. (2023). *Modelos de regressao beta para dados de escala*.
Master's dissertation, Universidade Federal do Parana, Curitiba. URI:
<https://hdl.handle.net/1884/86624>.

Ferrari, S. L. P., and Cribari-Neto, F. (2004). Beta regression for
modelling rates and proportions. *Journal of Applied Statistics*,
**31**(7), 799–815.
[doi:10.1080/0266476042000214501](https://doi.org/10.1080/0266476042000214501)

## See also

Useful links:

- <https://evandeilton.github.io/betaregscale/>

- <https://github.com/evandeilton/betaregscale>

- Report bugs at <https://github.com/evandeilton/betaregscale/issues>

## Author

**Maintainer**: José Evandeilton Lopes <evandeilton@gmail.com>
([ORCID](https://orcid.org/0009-0007-5887-4084))

Authors:

- Wagner Hugo Bonat ([ORCID](https://orcid.org/0000-0002-0349-7054))
