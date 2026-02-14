# betaregscale 2.1.0

## New features

* `betaregscale_simulate()` and `betaregscale_simulate_z()` gain a
  `delta` argument (default `NULL`) that forces all simulated
  observations to a specific censoring type: 0 (exact), 1 (left),
  2 (right), or 3 (interval). This enables targeted Monte Carlo
  studies with controlled censoring structures.
* `check_response()` gains a `delta` argument that accepts a vector
  of pre-specified censoring indicators, overriding the automatic
  boundary-based classification.

## Deprecations

* The `type` parameter (`"m"`, `"l"`, `"r"`) is deprecated across all
 functions: `betaregscale()`, `betaregscale_fit()`, `betaregscale_fit_z()`,
  `betaregscale_loglik()`, `betaregscale_loglik_z()`,
  `betaregscale_simulate()`, `betaregscale_simulate_z()`,
  `check_response()`, and `bs_prepare()`. Use `bs_prepare()` to
  control interval geometry instead. The parameter still works but
  emits a deprecation warning when passed explicitly.

# betaregscale 2.0.1

## New features

* **`bs_prepare()` data preprocessing**: new analyst-facing function that
  validates, classifies censoring, and rescales raw data before model fitting.
  Supports four flexible input modes: score-only, score + explicit delta,
  interval endpoints with NA patterns, and analyst-supplied left/right bounds.
  Prepared data is automatically detected by `betaregscale()`.
* Internal helper `.extract_response()` enables transparent detection of
  `bs_prepare()`-processed data across all fitting, log-likelihood, and
  starting-value functions.
* `censoring_summary()` now also accepts data frames from `bs_prepare()`.
* New vignette section documenting all four data preparation modes.

## Bug fixes

* Fixed potential row-indexing bug when `bs_prepare()` receives a subset
  data frame with non-sequential row names. Output now always has
  sequential row names (`1:n`).

# betaregscale 2.0.0

## Breaking changes

* Removed dependency on `bbmle`. All model fitting now uses `stats::optim()`
  directly with analytical gradients via the C++ backend.
* The `betaregscale_bbmle()` function has been removed.
* The `cumulative` parameter has been replaced by the `delta` indicator
  vector, which supports mixed censoring types within the same dataset.
* Parameter `dados` renamed to `data` across all functions.
* Simulation functions renamed: `betaregscale_simula_dados()` is now
  `betaregscale_simulate()`, and `betaregscale_simula_dados_z()` is now
  `betaregscale_simulate_z()`.

## New features

* **Mixed censoring support**: the complete likelihood (Eq. 2.24) now
  handles four censoring types simultaneously: exact ($\delta=0$),
  left-censored ($\delta=1$), right-censored ($\delta=2$), and
  interval-censored ($\delta=3$).
* **C++ backend rewrite**: log-likelihood and analytical gradient functions
  rewritten in C++ (RcppArmadillo) for numerically stable, high-performance
  evaluation.
* **betareg-style S3 interface**: `coef()` and `vcov()` now accept
  `model = c("full", "mean", "precision")` argument.
* New S3 methods: `nobs()`, `formula()`, `model.matrix()`, `confint()`,
  and `plot()`.
* `confint()` provides Wald confidence intervals based on the asymptotic
  normal approximation (z-test, not t-test).
* `plot()` method with six diagnostic panels (residuals vs indices,
  Cook's distance, residuals vs linear predictor, residuals vs fitted,
  half-normal envelope, predicted vs observed) and both base R and
  ggplot2 backends.
* `censoring_summary()` function for visual and tabular summaries of the
  censoring structure, with both base R and ggplot2 backends.
* `predict()` expanded with five types: `"response"`, `"link"`,
  `"precision"`, `"variance"`, and `"quantile"`. Supports `newdata` for
  both fixed and variable dispersion models.
* `residuals()` supports five types: `"response"`, `"pearson"`, `"rqr"`
  (randomized quantile residuals), `"weighted"`, and `"sweighted"`.
* `summary()` output now shows separate coefficient tables for mean and
  precision submodels with Wald z-tests.

## Bug fixes

* Fixed Pearson residual computation to correctly dispatch by
  reparameterization type (repar 1 vs repar 2).
* Fixed `predict()` with `newdata` for variable-dispersion models.
* Fixed p-values to use `pnorm()` (standard normal) instead of `pt()`
  (Student-t), consistent with Wald inference theory (Eq. 2.34--2.35).

# betaregscale 1.1.1

* Initial public release with `bbmle`-based fitting.
* Support for fixed and variable dispersion models.
* Basic S3 methods: `coef`, `vcov`, `fitted`, `residuals`, `summary`,
  `print`.
