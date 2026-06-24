# betaregscale 2.7.1

This is a maintenance release focused on correctness, numerical robustness, and
performance. It carries the fixes from a deep audit of the C++ backends and the R
interface. No user-facing API changes.

## Bug fixes

* Fixed the inverse link for the precision (dispersion) submodel, which could return a
  negative shape near the origin; the dispersion parameter is now always positive.
* Corrected deviance residuals to use the proper saturated-model log-likelihood (affecting
  `residuals()` and the diagnostic plots for both `brs` and `brsmm`). Negative discrepancies
  are now handled via `sqrt(abs(...))` rather than being silently truncated to zero.
* `predict.brs()` now detects a variable-dispersion model from the model's term labels
  instead of the number of parameters, fixing incorrect predictions for some fits.
* The `sqrt` link for the precision submodel now clamps the linear predictor to be
  non-negative, preserving the correct gradient sign during optimization.
* Mixed-effects mode finding (`brsmm`) is substantially more robust: the Newton-Raphson step
  is accepted only when it improves the objective, the Hessian regularization is now
  proportional to the smallest eigenvalue of `-H`, and non-finite proposed steps are rejected
  before evaluation, preventing `NaN`/`Inf` crashes.
* The adaptive Gauss-Hermite (AGHQ) integration grid now uses 64-bit indexing to prevent an
  integer-overflow crash for four or more random-effect dimensions, and the quasi-Monte Carlo
  prime table was expanded (20 to 50 primes) to keep Halton sequences uncorrelated in higher
  dimensions.
* Halton quantiles are clamped to `(1e-9, 1 - 1e-9)` to avoid infinite values at the boundary.

## Improvements

* The optimizer now emits a warning when it fails to converge in `brs()` and `brsmm()`.
* Pseudo-R2 reporting adds a caution note when more than 50% of observations are censored.
* Link evaluation no longer constructs `stats::make.link()` objects on every call:
  `apply_inv_link()` uses direct closed-form formulas, and a new `apply_link()` provides the
  forward transform, used across the fitting, methods, and plotting code.
* `brs_check()` is fully vectorized (no R-level loop) and `brs_prep()` uses `vapply()` with a
  vectorized `NA` guard.
* Starting values: `compute_start()` uses a moment-based estimate for the precision intercept
  (avoiding a second GLM fit), and `brsmm()` initializes the random-effect standard deviation
  from the between-group variance.
* `brs_coef()` now formally signals its deprecation via `.Deprecated("brs_est")`.
* `vcov.brs()` warns when it cannot compute a finite covariance matrix (e.g., when `MASS` is
  unavailable for the generalized-inverse fallback).

## Internal

* Extracted a shared C++ header (`src/brs_common.h`) to remove duplicated code between the
  Armadillo and Eigen backends, and documented the numerical tolerance constants.
* The Eigen backend reuses a single pre-allocated workspace vector when forming the numerical
  Hessian.
* `src/Makevars` and `src/Makevars.win` no longer define `-DARMA_NO_DEBUG`, enabling Armadillo
  bounds checking during development.
* The internal `.brsmm_loglik_eigen` entry point is exported with a leading dot to keep it out
  of the public namespace.

---

# betaregscale 2.7.0

## New features

* `autoplot.brsmm()` and `autoplot.brs()` gain three new arguments:
  * `theme`: accepts any ggplot2 theme object or function (default `ggplot2::theme_minimal()`),
    replacing the hardcoded theme in all 8 internal `brsmm` and 4 internal `brs` plot helpers.
  * `title`, `xlab`, `ylab`: override individual plot labels via `ggplot2::labs()`.
  * `type = "all"`: renders every available panel in a single `gridExtra::grid.arrange()` grid.
  * `ncol`: controls the number of columns when `type = "all"`.
  * `...`: passes additional named arguments to `ggplot2::theme()` on top of the base theme.
* `autoplot.brsmm()` adds `type = "shrinkage"`: scatter of Laplace posterior modes
  versus naïve per-group logit-mean deviations, with identity line and loess smoother.
* `autoplot.brsmm()` updates `type = "ranef_caterpillar"`: error bars now show ±1.96 × Model SD
  (marginal standard deviation from the `D` covariance matrix).
* `plot.brsmm()` adds two new base R diagnostic panels:
  * `which = 7`: Q-Q normal plot of random-effect posterior modes.
  * `which = 8`: dotchart caterpillar of posterior modes (ordered by value).
* `brsmm_re_study()` now returns `$icc`: intraclass correlation coefficient on the
  latent logistic scale (`σ²_b / (σ²_b + π²/3)`).
* `print.brsmm_re_study()` now displays a VarCorr-style table (Std.Dev., Corr) and
  the ICC alongside the existing shrinkage and normality diagnostics.

## Improvements

* Coefficient display: `print.summary.brsmm()`, `print.brsmm()`, `print.summary.brs()`,
  and `print.brs()` now apply cosmetic name cleaners that strip internal prefixes such as
  `(phi)_` from precision coefficients and convert Cholesky-factor internal names
  (e.g., `(re_chol_logsd)_X|g`) to readable `logSD.X|g` / `cov.X:Y|g` forms.
* Calibration plots (`autoplot.brs` and `autoplot.brsmm`) now map `linewidth = n`
  for `geom_line()` instead of the deprecated `size` aesthetic, eliminating the
  ggplot2 ≥ 3.4.0 deprecation warning.

## Testing

* Added `tests/testthat/test-re-and-autoplot-improvements.R` with 15 new tests covering:
  `.pretty_phi_names()`, `.pretty_re_names()`, `brsmm_re_study()` ICC and VarCorr output,
  `plot.brsmm()` panels 7 and 8, `autoplot.brsmm()` and `autoplot.brs()` title/xlab/ylab,
  theme arg (object and function), `type = "all"`, `...` forwarding.

---

# betaregscale 2.6.9


## CRAN resubmission

### Documentation and formatting fixes
* Added missing `\value`, `\seealso`, and `\examples{\donttest{...}}` tags to multiple S3 method documentation files (`print.summary`, `residuals`, `summary`, `vcov`, `ranef`) to ensure full CRAN policy compliance.
* Translated remaining Portuguese text into English in the mixed-effects vignette (`vignettes/brs-mm.Rmd`).
* Corrected `ranef()` usage in vignettes to correctly call the generic function.
* Fixed mathematical formulas rendering in `README.md` to be fully compatible with GitHub Markdown, and updated `pkgdown` site build configuration to load `betaregscale` appropriately during vignette setups.
* Minor mathematical formatting and typographical fixes (e.g., en-dashes for page ranges) in `README.md` references.# betaregscale 2.6.8

## New features
* Completed S3 method standardization for `brsmm` (mixed-effects) objects to mirror the interface of `brs` (fixed-effects) objects:
    * Added missing extractors: `formula()`, `model.matrix()`, and `confint()`.
    * Upgraded `residuals()` to support conditional `"deviance"`, `"rqr"` (randomized quantile residuals), `"weighted"`, and `"sweighted"` options.
    * Upgraded `predict()` to support conditional `type = "quantile"` evaluations directly.
    * Added `ranef()` generic and `ranef.brsmm()` method to extract random-effect modes.
* Modified package helper functions `brs_gof()` and `brs_est()` to compute GOF properties and estimates directly from both `brs` and `brsmm` objects respectively.

## Improvements
* Standardized `print.brsmm()` to explicitly display mean, precision, and random-effect coefficient blocks side-by-side, mirroring the verbose visual style of `print.brs()`.

## Documentation
* **Complete `@examples` audit**: added runnable `\donttest{}` examples to all ~30 previously undocumented exported functions, including all S3 methods for `brs` and `brsmm` objects (`coef()`, `vcov()`, `logLik()`, `AIC()`, `BIC()`, `nobs()`, `formula()`, `model.matrix()`, `fitted()`, `residuals()`, `confint()`, `predict()`, `print()`, `summary()`, `ranef()`, `anova()`, `plot()`, `autoplot()`).
* **Removed all `set.seed()` calls from examples** across 15+ files (`fit.R`, `brsmm.R`, `bootstrap.R`, `cv.R`, `marginaleffects.R`, `scoreprob.R`, `table.R`, `simulate.R`, `prepare.R`, `autoplot.R`, `autoplot-brsmm.R`, `loglik.R`). All examples now use deterministic toy datasets.
* **No `\dontrun{}` anywhere**: all examples are either direct or wrapped in `\donttest{}` as appropriate.
* **Ferrari & Cribari-Neto (2004) DOI** (`10.1080/0266476042000214501`) added to every occurrence of that reference across `betaregscale-package.R`, `brsmm.R`, `methods.R`, `anova-methods.R`, and `brsmm-random-effects-study.R`.
* **`@seealso` cross-links** added to all S3 method documentation blocks for both `brs` and `brsmm` objects.
* **`brs_coef()`** documentation updated with deprecation notice, `@description`, `@return`, and `@seealso`.
* **`brs_hessian()`** documentation improved: added `@param object`, `@seealso`, and a deterministic example.
* **`print.brsmm_re_study()`** now has a complete roxygen2 block including `@description`, `@param`, `@return`, `@method`, `@seealso`, and `@examples`.
* **`ranef()` generic** now includes `@param`, `@return`, `@seealso`, and `@examples`.
* All `autoplot.*` examples updated to use `ggplot2::autoplot()` (explicit namespace) for reliability in check environments.

---

# betaregscale 2.6.7
## CRAN resubmission (Konstanze Lauseker review, 20 Feb 2026)

### Bug fixes and CRAN policy compliance

* Added `\value` documentation to `print.brs()` and `print.summary.brs()` methods.
* Replaced `\dontrun{}` with `\donttest{}` in `brs_gof()` example and created complete executable example.
* Removed `.GlobalEnv` modification from `brs_bootstrap()` (CRAN policy violation).
* Removed `set.seed()` calls from exported functions: `brs_bootstrap()`, `brs_marginaleffects()`, `brsmm()`, and `brs_cv()`. Users must now call `set.seed()` externally before these functions for reproducibility.
* Removed `seed` parameter from all four functions listed above. Documentation updated with recommended usage pattern.

---

# betaregscale 2.6.6

## CRAN resubmission (Uwe Ligges review, 18 Feb 2026)

* DESCRIPTION: function names in Title/Description now use parentheses (e.g. `logLik()`, `coef()`) per CRAN policy.
* URL updated to GitHub repository (https://github.com/evandeilton/betaregscale) while pkgdown site is deployed.

---

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
