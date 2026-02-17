# betaregscale 2.6.5

## New features

* Extended `brs_bootstrap()` with `ci_type = "bca"` (bias-corrected and accelerated intervals), plus Monte Carlo diagnostics for interval endpoints (`mcse_lower`, `mcse_upper`).
* Added Wald interval columns (`wald_lower`, `wald_upper`) to bootstrap output for direct asymptotic vs resampling comparison.
* Added `autoplot.brs_bootstrap()` support to visually compare bootstrap and Wald intervals in `type = "ci_forest"`.
* Added `autoplot.brs_marginaleffects()` with three views: `forest`, `magnitude`, and `dist`.

## Improvements

* Improved robustness and efficiency in `brs_marginaleffects()`:
  - central-difference AME approximation for numeric covariates,
  - scale-adaptive perturbation step,
  - one-time simulation draw generation reused across variables,
  - optional storage of AME draws via `keep_draws = TRUE`.
* Refined `brs_cens()` output to include richer summary fields (`percentage`, `severity`, `interpretation`) and optional domain-agnostic interpretation messages via `inform = TRUE`.
* Updated README and vignettes with examples for BCa bootstrap intervals, bootstrap visual diagnostics, and enhanced marginal-effects visualization workflow.

---

# betaregscale 2.6.4

## New features

* Extended `brsmm()` to support multivariate random effects in the mean predictor, including random intercept + random slope specifications such as `random = ~ 1 + x | group`.
* Added multivariate Laplace approximation in the Eigen C++ backend for group-specific latent vectors and covariance matrix handling via packed lower-Cholesky parameterization.
* Added `brsmm_group_modes_eigen()` to compute posterior modes of group random effects for general random-effects dimension.
* Added generic model-comparison methods `anova.brs()` and `anova.brsmm()` for likelihood-ratio workflow across `brs` and `brsmm` candidates.
* Added `brsmm_re_study()` and `print.brsmm_re_study()` for numeric random-effects diagnostics (covariance/correlation, shrinkage, normality checks).

## Improvements

* Updated `predict.brsmm()`, `vcov.brsmm()`, and `print.brsmm()` to support both scalar (`q_b=1`) and vector (`q_b>1`) random-effects structures.
* Expanded mixed-model test coverage with integration tests for random intercept + slope fits, covariance extraction (`D`), `ranef`, random-effects studies, and prediction behavior.
* Updated `README.md` and vignettes with explicit multivariate mixed-model mathematics, Laplace formula in matrix form, and end-to-end model-selection examples.
* Documentation references now use DOI-based validated links only (`https://doi.org/...`) to keep CRAN URL checks robust.

---

# betaregscale 2.6.3

## Improvements

* Revised and expanded all core vignettes (`brs-intro`, `brs-analyst-tools`, `brs-mm`) with stronger mathematical exposition, explicit likelihood pieces by censoring type, and clearer inferential interpretation for analysts.
* Updated vignettes and README to prioritize the package's most important analyst-facing functions: `brs_bootstrap()`, `brs_marginaleffects()`, `brs_predict_scoreprob()`, `brs_cv()`, and `brs_table()`.
* Standardized vignette outputs with cleaner tabular presentation using `knitr::kable(..., digits = 4)` for better readability and reporting consistency.
* Added and revised bibliographic references with validated DOI metadata and dual online source verification links in vignettes/README.
* Re-rendered all vignettes and rebuilt documentation website (`pkgdown`) to keep articles and reference pages synchronized with the current API.

---

# betaregscale 2.6.2

## Improvements

* Improved numerical stability in `brsmm()` by refining the optimization control and starting values.
* Updated `simulate()` method to better handle edge cases in random effects simulation.
* Enhanced `methods.R` for better compatibility with downstream packages.

---

# betaregscale 2.6.1

## Bug fixes

* Renamed vignettes to avoid naming collisions with the package name, which caused `pkgdown` site build failures.
* Updated `_pkgdown.yml` to reflect new vignette names.

---

# betaregscale 2.6.0

## New features

* Added `brsmm()` for mixed-effects beta interval regression with
  Gaussian random intercepts (`random = ~ 1 | group`) using
  Laplace-approximated marginal likelihood.
* Added C++ mixed-model likelihood core:
  `.brsmm_loglik_laplace_cpp()` and `.brsmm_group_modes_cpp()`.
* Added a first S3 interface for `brsmm` objects:
  `print`, `summary`, `coef`, `vcov`, `logLik`, `AIC`, `BIC`,
  `nobs`, `fitted`, `predict`, and `residuals`.

## Improvements

* Added `test-brsmm.R` with mixed-model fitting and prediction tests.
* Corrected author name spelling in package metadata/documentation:
  **José Evandeilton Lopes**.

---

# betaregscale 2.5.0

## New features

* Added `brs_table()` to compare one or more fitted `brs` models in a
  single table with `logLik`, `AIC`, `BIC`, pseudo-R2, and censoring
  composition.
* Added `brs_marginaleffects()` for average marginal effects in the
  mean or precision submodel, with optional simulation-based
  uncertainty intervals.
* Added `autoplot.brs()` with `ggplot2` diagnostics for
  `type = "calibration"`, `type = "score_dist"`, `type = "cdf"`, and
  `type = "residuals_by_delta"`.
* Added `brs_predict_scoreprob()` to obtain predicted probabilities on
  the original integer score scale.
* Added `brs_cv()` for repeated k-fold cross-validation of `brs`
  models with fold-level predictive metrics (`log_score`, `rmse_yt`,
  and `mae_yt`).

## Improvements

* Updated package reference organization (`pkgdown`) to expose the new
  analyst-oriented tools.
* Updated `README.md` and vignette content with examples for model
  comparison, marginal effects, and score-probability predictions.

---

# betaregscale 2.4.0

## Breaking changes

* `brs_sim_var()` is no longer exported. Variable-dispersion simulation is now done through `brs_sim()` using a two-part formula (for example, `~ x1 + x2 | z1 + z2`).
* `brs_loglik()` and `brs_loglik_var()` are now internal helpers and are no longer part of the user-facing API.

## New features

* `brs_sim()` is now the single simulation entry point for both fixed- and variable-dispersion models, with formula semantics aligned to `brs()`.

## Improvements

* Release documentation was updated to reflect the consolidated simulation API and current exported function set.
* `brs_prep()` consistency warnings are emitted once per call on final prepared output, improving test stability and warning capture behavior.

---

# betaregscale 2.3.0

## Breaking changes

* **API Overhaul**: All exported functions have been renamed to use the compact `brs_` prefix for consistency and ease of typing.
    * `betaregscale()` -> `brs()`
    * `betaregscale_fit()` -> `brs_fit_fixed()`
    * `betaregscale_fit_z()` -> `brs_fit_var()`
    * `betaregscale_loglik()` -> `brs_loglik()`
    * `betaregscale_loglik_z()` -> `brs_loglik_var()`
    * `betaregscale_simulate()` -> `brs_sim()`
    * `betaregscale_simulate_z()` -> `brs_sim_var()`
    * `prepare_data()` -> `brs_prep()`
    * `check_response()` -> `brs_check()`
    * `censoring_summary()` -> `brs_cens()`
    * `beta_reparam()` -> `brs_repar()`
    * `gof()` -> `brs_gof()`
    * `est()` -> `brs_est()`
    * `hessian_matrix()` -> `brs_hessian()`
    * `betaregscale_coef()` -> `brs_coef()`

* **Class Renaming**: The S3 class `betaregscale` has been renamed to `brs`. All associated S3 methods have been updated accordingly (e.g., `summary.brs`, `plot.brs`).

---

# betaregscale 2.2.0

## Breaking changes

* **`type` argument removed**: The deprecated `type` argument has been
  completely removed from all functions: `check_response()`,
  `prepare_data()`, `betaregscale()`, `betaregscale_fit()`,
  `betaregscale_fit_z()`, `betaregscale_loglik()`,
  `betaregscale_loglik_z()`, `betaregscale_simulate()`,
  `betaregscale_simulate_z()`, and internal helpers `compute_start()`,
  `.extract_response()`, `.build_simulated_response()`, and
  `.compute_endpoints()`.  The midpoint interval geometry
  (`type = "m"`) is now the only option and is hardcoded internally.
  Users who previously relied on `type = "l"` or `type = "r"` should
  use `prepare_data()` to supply custom left/right endpoints instead.

* **Renamed `bs_prepare()` to `prepare_data()`**: The data preparation
  function has been renamed to `prepare_data()` to be more descriptive
  and consistent with the package's verb-based API. The returned data
  frame now carries the `is_prepared` attribute instead of `bs_prepared`.

---

# betaregscale 2.1.1

## New features

* **`delta` argument in simulation functions**:
  `betaregscale_simulate()` and `betaregscale_simulate_z()` gain a
  `delta` argument (default `NULL`) that forces all simulated
  observations to a specific censoring type: 0 (exact), 1 (left),
  2 (right), or 3 (interval).  This enables targeted Monte Carlo
  studies where the analyst controls the censoring structure.

  When `delta` is non-NULL, the actual simulated values
  (`y_raw = rbeta(n, a, b)`) are preserved on the scale grid, and
  the forced censoring indicator is passed to `check_response()` as
  a vector.  This ensures that each observation retains its
  covariate-driven variation with observation-specific endpoints.

  The returned data frame carries `attr(, "bs_prepared") = TRUE` so
  that `betaregscale()`, `betaregscale_loglik()`, and all fitting
  functions use the pre-computed `left`, `right`, `yt`, and `delta`
  columns directly, bypassing the automatic boundary classification.
  Without this attribute, the fitting pipeline would re-classify the
  response from the `y` column alone, which would ignore the forced
  delta.

* **`delta` argument in `check_response()`**: accepts an integer
  vector of pre-specified censoring indicators, overriding the
  automatic boundary-based classification on a per-observation
  basis.  The endpoint formulas adapt to non-boundary observations:

  | delta | condition | left (l_i)      | right (u_i)     |
  |-------|-----------|-----------------|-----------------|
  | 0     | any       | y / K           | y / K           |
  | 1     | y = 0     | eps             | lim / K         |
  | 1     | y != 0    | eps             | (y + lim) / K   |
  | 2     | y = K     | (K - lim) / K   | 1 - eps         |
  | 2     | y != K    | (y - lim) / K   | 1 - eps         |
  | 3     | type "m"  | (y - lim) / K   | (y + lim) / K   |

  The distinction between boundary and non-boundary observations is
  essential: when delta = 1 is forced on a non-zero y, the upper
  bound uses the actual y value ((y + lim)/K) rather than the fixed
  boundary formula (lim/K).  This preserves the information content
  of each observation.

* **Observation-specific endpoints in `bs_prepare()`**: the internal
  `.compute_endpoints()` helper now uses the same adaptive formulas
  as `check_response()` for analyst-forced left/right censoring on
  non-boundary scores.  Previously, delta = 1 always produced
  `right = lim/K` and delta = 2 always produced
  `left = (K - lim)/K`, regardless of the actual y value.

## Bug fixes

* **Simulation with forced `delta = 1` or `delta = 2`**: the
  internal `.build_simulated_response()` helper previously replaced
  all y values with boundary values (`y_grid = rep(0, n)` for
  delta = 1, `y_grid = rep(ncuts, n)` for delta = 2).  This
  produced degenerate data where every observation had identical
  endpoints (e.g., all `left = 0.995, right = 0.99999` for
  delta = 2), destroying all covariate-driven variation and making
  regression fitting impossible.

  The fix preserves the actual simulated grid values
  (`y_grid = round(y_raw * ncuts)`) and passes a forced delta
  vector to `check_response()`, which computes observation-specific
  endpoints using the actual y values.

* **Missing `"bs_prepared"` attribute on simulation output**: when
  `delta` was forced, the simulation functions did not mark the
  output with `attr(, "bs_prepared") = TRUE`.  As a result,
  `betaregscale()` would re-classify the response via
  `check_response()`, silently overwriting the forced delta with
  automatic boundary rules.  The attribute is now set correctly.

## Deprecations

* The `type` parameter (`"m"`, `"l"`, `"r"`) is deprecated across all
  functions: `betaregscale()`, `betaregscale_fit()`,
  `betaregscale_fit_z()`, `betaregscale_loglik()`,
  `betaregscale_loglik_z()`, `betaregscale_simulate()`,
  `betaregscale_simulate_z()`, `check_response()`, and
  `prepare_data()`. Use `prepare_data()` to control interval geometry
  instead. The parameter still works but emits a deprecation warning
  when passed explicitly.

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
