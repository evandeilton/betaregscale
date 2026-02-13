# Changelog

## betaregscale 2.0.1

### New features

- **[`bs_prepare()`](https://evandeilton.github.io/betaregscale/reference/bs_prepare.md)
  data preprocessing**: new analyst-facing function that validates,
  classifies censoring, and rescales raw data before model fitting.
  Supports four flexible input modes: score-only, score + explicit
  delta, interval endpoints with NA patterns, and analyst-supplied
  left/right bounds. Prepared data is automatically detected by
  [`betaregscale()`](https://evandeilton.github.io/betaregscale/reference/betaregscale.md).
- Internal helper `.extract_response()` enables transparent detection of
  [`bs_prepare()`](https://evandeilton.github.io/betaregscale/reference/bs_prepare.md)-processed
  data across all fitting, log-likelihood, and starting-value functions.
- [`censoring_summary()`](https://evandeilton.github.io/betaregscale/reference/censoring_summary.md)
  now also accepts data frames from
  [`bs_prepare()`](https://evandeilton.github.io/betaregscale/reference/bs_prepare.md).
- New vignette section documenting all four data preparation modes.

### Bug fixes

- Fixed potential row-indexing bug when
  [`bs_prepare()`](https://evandeilton.github.io/betaregscale/reference/bs_prepare.md)
  receives a subset data frame with non-sequential row names. Output now
  always has sequential row names (`1:n`).

## betaregscale 2.0.0

### Breaking changes

- Removed dependency on `bbmle`. All model fitting now uses
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html) directly with
  analytical gradients via the C++ backend.
- The `betaregscale_bbmle()` function has been removed.
- The `cumulative` parameter has been replaced by the `delta` indicator
  vector, which supports mixed censoring types within the same dataset.
- Parameter `dados` renamed to `data` across all functions.
- Simulation functions renamed: `betaregscale_simula_dados()` is now
  [`betaregscale_simulate()`](https://evandeilton.github.io/betaregscale/reference/betaregscale_simulate.md),
  and `betaregscale_simula_dados_z()` is now
  [`betaregscale_simulate_z()`](https://evandeilton.github.io/betaregscale/reference/betaregscale_simulate_z.md).

### New features

- **Mixed censoring support**: the complete likelihood (Eq. 2.24) now
  handles four censoring types simultaneously: exact ($\delta = 0$),
  left-censored ($\delta = 1$), right-censored ($\delta = 2$), and
  interval-censored ($\delta = 3$).
- **C++ backend rewrite**: log-likelihood and analytical gradient
  functions rewritten in C++ (RcppArmadillo) for numerically stable,
  high-performance evaluation.
- **betareg-style S3 interface**:
  [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) now accept
  `model = c("full", "mean", "precision")` argument.
- New S3 methods: [`nobs()`](https://rdrr.io/r/stats/nobs.html),
  [`formula()`](https://rdrr.io/r/stats/formula.html),
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html), and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).
- [`confint()`](https://rdrr.io/r/stats/confint.html) provides Wald
  confidence intervals based on the asymptotic normal approximation
  (z-test, not t-test).
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method with
  six diagnostic panels (residuals vs indices, Cook’s distance,
  residuals vs linear predictor, residuals vs fitted, half-normal
  envelope, predicted vs observed) and both base R and ggplot2 backends.
- [`censoring_summary()`](https://evandeilton.github.io/betaregscale/reference/censoring_summary.md)
  function for visual and tabular summaries of the censoring structure,
  with both base R and ggplot2 backends.
- [`predict()`](https://rdrr.io/r/stats/predict.html) expanded with five
  types: `"response"`, `"link"`, `"precision"`, `"variance"`, and
  `"quantile"`. Supports `newdata` for both fixed and variable
  dispersion models.
- [`residuals()`](https://rdrr.io/r/stats/residuals.html) supports five
  types: `"response"`, `"pearson"`, `"rqr"` (randomized quantile
  residuals), `"weighted"`, and `"sweighted"`.
- [`summary()`](https://rdrr.io/r/base/summary.html) output now shows
  separate coefficient tables for mean and precision submodels with Wald
  z-tests.

### Bug fixes

- Fixed Pearson residual computation to correctly dispatch by
  reparameterization type (repar 1 vs repar 2).
- Fixed [`predict()`](https://rdrr.io/r/stats/predict.html) with
  `newdata` for variable-dispersion models.
- Fixed p-values to use [`pnorm()`](https://rdrr.io/r/stats/Normal.html)
  (standard normal) instead of
  [`pt()`](https://rdrr.io/r/stats/TDist.html) (Student-t), consistent
  with Wald inference theory (Eq. 2.34–2.35).

## betaregscale 1.1.1

- Initial public release with `bbmle`-based fitting.
- Support for fixed and variable dispersion models.
- Basic S3 methods: `coef`, `vcov`, `fitted`, `residuals`, `summary`,
  `print`.
