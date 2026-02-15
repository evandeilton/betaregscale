# Changelog

## betaregscale 2.1.1

### New features

- **`delta` argument in simulation functions**:
  [`betaregscale_simulate()`](https://evandeilton.github.io/betaregscale/reference/betaregscale_simulate.md)
  and
  [`betaregscale_simulate_z()`](https://evandeilton.github.io/betaregscale/reference/betaregscale_simulate_z.md)
  gain a `delta` argument (default `NULL`) that forces all simulated
  observations to a specific censoring type: 0 (exact), 1 (left), 2
  (right), or 3 (interval). This enables targeted Monte Carlo studies
  where the analyst controls the censoring structure.

  When `delta` is non-NULL, the actual simulated values
  (`y_raw = rbeta(n, a, b)`) are preserved on the scale grid, and the
  forced censoring indicator is passed to
  [`check_response()`](https://evandeilton.github.io/betaregscale/reference/check_response.md)
  as a vector. This ensures that each observation retains its
  covariate-driven variation with observation-specific endpoints.

  The returned data frame carries `attr(, "bs_prepared") = TRUE` so that
  [`betaregscale()`](https://evandeilton.github.io/betaregscale/reference/betaregscale.md),
  [`betaregscale_loglik()`](https://evandeilton.github.io/betaregscale/reference/betaregscale_loglik.md),
  and all fitting functions use the pre-computed `left`, `right`, `yt`,
  and `delta` columns directly, bypassing the automatic boundary
  classification. Without this attribute, the fitting pipeline would
  re-classify the response from the `y` column alone, which would ignore
  the forced delta.

- **`delta` argument in
  [`check_response()`](https://evandeilton.github.io/betaregscale/reference/check_response.md)**:
  accepts an integer vector of pre-specified censoring indicators,
  overriding the automatic boundary-based classification on a
  per-observation basis. The endpoint formulas adapt to non-boundary
  observations:

  | delta | condition | left (l_i)    | right (u_i)   |
  |-------|-----------|---------------|---------------|
  | 0     | any       | y / K         | y / K         |
  | 1     | y = 0     | eps           | lim / K       |
  | 1     | y != 0    | eps           | (y + lim) / K |
  | 2     | y = K     | (K - lim) / K | 1 - eps       |
  | 2     | y != K    | (y - lim) / K | 1 - eps       |
  | 3     | type “m”  | (y - lim) / K | (y + lim) / K |

  The distinction between boundary and non-boundary observations is
  essential: when delta = 1 is forced on a non-zero y, the upper bound
  uses the actual y value ((y + lim)/K) rather than the fixed boundary
  formula (lim/K). This preserves the information content of each
  observation.

- **Observation-specific endpoints in
  [`bs_prepare()`](https://evandeilton.github.io/betaregscale/reference/bs_prepare.md)**:
  the internal `.compute_endpoints()` helper now uses the same adaptive
  formulas as
  [`check_response()`](https://evandeilton.github.io/betaregscale/reference/check_response.md)
  for analyst-forced left/right censoring on non-boundary scores.
  Previously, delta = 1 always produced `right = lim/K` and delta = 2
  always produced `left = (K - lim)/K`, regardless of the actual y
  value.

### Bug fixes

- **Simulation with forced `delta = 1` or `delta = 2`**: the internal
  `.build_simulated_response()` helper previously replaced all y values
  with boundary values (`y_grid = rep(0, n)` for delta = 1,
  `y_grid = rep(ncuts, n)` for delta = 2). This produced degenerate data
  where every observation had identical endpoints (e.g., all
  `left = 0.995, right = 0.99999` for delta = 2), destroying all
  covariate-driven variation and making regression fitting impossible.

  The fix preserves the actual simulated grid values
  (`y_grid = round(y_raw * ncuts)`) and passes a forced delta vector to
  [`check_response()`](https://evandeilton.github.io/betaregscale/reference/check_response.md),
  which computes observation-specific endpoints using the actual y
  values.

- **Missing `"bs_prepared"` attribute on simulation output**: when
  `delta` was forced, the simulation functions did not mark the output
  with `attr(, "bs_prepared") = TRUE`. As a result,
  [`betaregscale()`](https://evandeilton.github.io/betaregscale/reference/betaregscale.md)
  would re-classify the response via
  [`check_response()`](https://evandeilton.github.io/betaregscale/reference/check_response.md),
  silently overwriting the forced delta with automatic boundary rules.
  The attribute is now set correctly.

### Deprecations

- The `type` parameter (`"m"`, `"l"`, `"r"`) is deprecated across all
  functions:
  [`betaregscale()`](https://evandeilton.github.io/betaregscale/reference/betaregscale.md),
  [`betaregscale_fit()`](https://evandeilton.github.io/betaregscale/reference/betaregscale_fit.md),
  [`betaregscale_fit_z()`](https://evandeilton.github.io/betaregscale/reference/betaregscale_fit_z.md),
  [`betaregscale_loglik()`](https://evandeilton.github.io/betaregscale/reference/betaregscale_loglik.md),
  [`betaregscale_loglik_z()`](https://evandeilton.github.io/betaregscale/reference/betaregscale_loglik_z.md),
  [`betaregscale_simulate()`](https://evandeilton.github.io/betaregscale/reference/betaregscale_simulate.md),
  [`betaregscale_simulate_z()`](https://evandeilton.github.io/betaregscale/reference/betaregscale_simulate_z.md),
  [`check_response()`](https://evandeilton.github.io/betaregscale/reference/check_response.md),
  and
  [`bs_prepare()`](https://evandeilton.github.io/betaregscale/reference/bs_prepare.md).
  Use
  [`bs_prepare()`](https://evandeilton.github.io/betaregscale/reference/bs_prepare.md)
  to control interval geometry instead. The parameter still works but
  emits a deprecation warning when passed explicitly.

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
