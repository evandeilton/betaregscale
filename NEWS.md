# betaregscale 3.0.0

Major release after an audit of the R code and the compiled backend
(2026-09), checked against Lopes (2023). The major version reflects that some
calls accepted by 2.7.4 now stop and some results change. With the defaults
(`repar = 2`, `interval = "mid"`, `lim = 0.5`) and well-behaved data,
estimates agree with 2.7.4 to optimiser tolerance; the sections below give
every change and its reason.

## Breaking changes

* Links are checked against `repar`: `identity`, `inverse` and `1/mu^2` are
  rejected for positive parameters, and under `repar = 1` the default
  `link_phi` is now `"log"` (it was `"logit"`, which capped the precision
  at 1).
* Fits change on data with observations far in a tail (the `1e-15`
  probability floor of the likelihood is gone), and new starting values move
  some variable-dispersion fits to higher optima.
* Input accepted with a warning or silently now stops: `lim` outside
  `(0, 0.5]`, scores outside `0..ncuts`, `delta = 3` with `left == right`,
  `NA` in a prepared `delta`, rank-deficient designs.
* `anova()` refuses fits with different `interval`, `ncuts`, `lim` or `nobs`.
* `summary()` of a `brsmm` fit reports SD/Corr intervals instead of z-tests
  on `log(sd)`; the ICC of `brsmm_re_study()` uses the beta level-1 variance
  instead of `pi^2/3`.
* `vcov()` returns `NA` (it returned 0) for variances it cannot estimate. The
  default Hessian is the compiled one (`hessian_method = "cpp"`; standard
  errors agree with `numDeriv` to 1e-8).
* `brs_check()` treats values in `(0, 1)` as exact per observation, as
  `brs_prep()` does.
* QMC with two or more random effects is a different estimator
  (symmetric-root scaling of the nodes).

## Bootstrap, fit diagnostics and mixed-model inference

* `brs_bootstrap()` simulates only the response, at the fitted shapes, and
  refits the original formula on a copy of the data: factors, `log(x)`,
  `0 + x`, any response name and variable dispersion work (they ended in
  "Too few successful bootstrap replicates (0)"). Each row keeps its
  observation mechanism: exact values stay continuous (they were re-gridded
  to `delta = 3`), scores are re-coarsened on the fit's grid, analyst
  thresholds (`brs_prep()` Modes 2-4) are kept as fixed, non-informative
  thresholds and `delta` is re-drawn by the cell of their partition where the
  new value falls (conservative when the original design had more thresholds
  than a row records; an approximation for Mode 2 rows with a forced `delta`).
  Failed replicates are counted (`n_failed`, `fail_rate`) and printed; the
  "basic" MCSE limits were swapped; `"bca"` warns once per session that it is
  an approximation.
* `brs()`/`brsmm()` stop on rank-deficient design matrices, warn on near
  collinearity, and check the gradient, the Hessian and the likelihood
  clamps after `optim()` (`fit$diagnostics`, one-line warnings; the gradient
  check is the log-likelihood gain of the remaining Newton step, above 0.01).
  `vcov()` no longer uses a generalised inverse: variances it cannot estimate
  are `NA` (they were 0). `brsmm()` flags a variance component on the
  boundary. `logLik()` of a `brs()` fit is evaluated exactly at the returned
  estimate.
* `brsmm`: `summary()` reports SD/Corr with transformed Wald intervals and
  no z-test on `log(sd)`; `anova()` uses the chi-bar-square mixture
  `1/2 chi2(Df - 1) + 1/2 chi2(Df)` when a model adds one random-effect term
  (`1/2 chi2(0) + 1/2 chi2(1)` for `brs` against a random intercept) and says
  so in its heading; the ICC of `brsmm_re_study()` uses the beta level-1
  variance of `logit(Y)` instead of `pi^2/3` (`NA` with a warning when the
  clamp of the mean drives it, e.g. probit or cloglog links with a large
  random-effect SD); `control` is merged into `list(maxit = 2000)`.
* `summary()` leaves the RNG state untouched. `brs_check()` treats values in
  `(0, 1)` as exact per observation, as `brs_prep()` does, and warns on input
  mixing them with values `>= 1` (half-point scores: use `y * 2` and
  `ncuts * 2`). `brs_sim(delta = 1 or 2)` warns
  (informative censoring, no finite MLE) and keeps covariates of `0 + x`
  formulas; `brs_prep()` warns on rows covering the whole scale. Marginal
  effects honour `| 0 + z`; repar 2 variable-dispersion starts use the
  moment intercept with zero slopes. This changes some fits: in a validation
  battery 22 of 136 fits moved, all to an equal or higher log-likelihood with
  fewer iterations, and a `brsmm` variable-dispersion fit (200 groups x 25)
  that had "converged" to a wrong point (intercept 0.94, true 0.2) now gains
  1214 in log-likelihood (intercept 0.14).

## Interval direction (`interval`)

* New `interval = c("mid", "right", "left")` in `brs_check()`, `brs_prep()`,
  `brs()`, `brsmm()` and `brs_sim()` (the dissertation's `m`, `r`, `l`). Score
  `s` maps to `[s - lim, s + lim] / K` (`"mid"`, default, unchanged) or to
  `[s, s + 1] / (K + 1)` (`"right"`/`"left"`: `K + 1` equal cells). This is a
  package choice that differs from the dissertation (which divides by `K`):
  `"right"` and `"left"` give the same likelihood and differ only in
  `predict(type = "score")`, so chapter 4's opposite intercept biases vanish.
  `delta` comes from the score (`0 -> 1`, `K -> 2`, else `3`).

* `brs_sim()` generates scores by the same mechanism the likelihood assumes
  and always attaches `is_prepared`/`ncuts`/`lim`/`interval`; `brs()` and
  `brsmm()` reuse these attributes.

* `lim` outside `(0, 0.5]` is an error, `lim < 0.5` warns (partial
  coarsening). `anova()` refuses fits with different `interval`, `ncuts`,
  `lim` (only under `"mid"`) or `nobs`. New `predict()` types `"score"`
  (latent score) and `"expected_score"` (expected recorded score).

* `brs_check()` now stops on scores outside `0..K` (it warned and the
  observation contributed `-1e6` to the log-likelihood), and `brs_check()`/
  `brs_prep()` reject `delta = 3` with `left == right`. `brs_prep()` checks
  analyst bounds against the latent range of the direction and makes an
  analyst interval that reaches 0 or 1 left- or right-censored.

## Reparameterizations and links

* The three schemes of `brs_repar()` model different parameters, so the
  admissible links now depend on `repar` and incompatible combinations are
  rejected with an error that prints the table: `repar = 0` (shapes `p, q`,
  both on `(0, Inf)`): `link` and `link_phi` in `{log, sqrt}`; `repar = 1`
  (mean, precision): `link` in `{logit, probit, cauchit, cloglog}`,
  `link_phi` in `{log, sqrt}`; `repar = 2` (mean, dispersion): both in
  `{logit, probit, cauchit, cloglog}`. `identity`, `inverse` and `1/mu^2`
  are no longer accepted for positive parameters (their inverse does not map
  the real line onto `(0, Inf)`); with `sqrt` a warning is issued after the
  fit when a fitted linear predictor is `<= 0` (flat inverse link). Link
  names are matched exactly (no partial matching).

* `link` and `link_phi` default to `NULL` in `brs()`, `brs_fit_fixed()`,
  `brs_fit_var()`, `brsmm()` and `brs_sim()`, and resolve by `repar`:
  `0 -> log/log`, `1 -> logit/log`, `2 -> logit/logit`. The old default
  `link_phi = "logit"` under `repar = 1` squashed the precision into
  `(0, 1)`: with true precision 20 it estimated 0.999999 (log-likelihood
  -2483 against -1885 with the `log` link).

* Under `repar = 0` the fitted object keeps the shape `p` in `hatmu` /
  `fitted_mu` (and `brs_repar(mu = )` takes `p > 0`), but every user-facing
  mean now is `E[Y] = a / (a + b)`: `fitted()`, `predict(type =
  "response")`, response/Pearson residuals, marginal effects and the
  calibration plots. `predict()` used to return `p` (0.949 where the mean
  was 0.493). Weighted/sweighted residuals are computed from the shapes
  (`digamma(a) - digamma(b)`, `trigamma(a) + trigamma(b)`), the pseudo
  R-squared compares `E[Y]` and `y` on the logit scale under `repar = 0`,
  and starting values for the shapes come from the method of moments.
  `brs_predict_scoreprob(newdata = )`, `brs_cv()` and the CDF plot no
  longer feed `predict(type = "response")` back into `brs_repar()`.
  Identical results under `repar = 1, 2`.

* `brsmm()` now runs the same input validation as `brs()`: `repar` outside
  `0:2` (which used to fall into the C++ `default:` branch silently),
  `ncuts`, `lim` and the links.

* `ncuts` and `lim` default to `NULL` in `brs()` and `brsmm()` and are taken
  from the attributes that `brs_prep()` (and now `brs_sim()`) store on the
  data; an explicit value that differs from the stored one is ignored with
  a warning. Refitting prepared data with another `ncuts` made
  `brs_predict_scoreprob()` rows sum to 0.05 instead of 1. `brs_sim()`
  always attaches `is_prepared`, `ncuts` and `lim` (it did so only when
  `delta` was given).

* The compiled likelihood clamped the inverse-linked first parameter to
  `(1e-5, 1 - 1e-5)` for every `repar`, so under `repar = 0` the truth and
  the "MLE" were both evaluated with the shape `p` capped at `0.99999`.
  `clamp_mu_by_repar()` in `src/brs_common.h` now uses `[1e-5, 1e8]` for
  the shape and `[1e-5, 1 - 1e-5]` for the mean, and the R side (`hatmu`,
  `hatphi`, `predict()`) applies the same clamps.

## Row alignment and input validation

* `.extract_response()` and `.brsmm_row_index()` treated numeric row names as
  row positions. After `data[-10, ]`, a permutation, or any subset that keeps
  the original row names, `left`/`right`/`delta` and the grouping variable
  were taken from the wrong rows, silently corrupting `brs()` on subsetted
  data, `brs_cv()`, the BCa jackknife in `brs_bootstrap()` and `brsmm()`.
  Rows are now always mapped by `match(rownames(mf), rownames(data))`.

* The compiled `brsmm()` likelihood did no dimension checks. With a `NA` in a
  random-slope variable, `model.matrix()` dropped a row and the group builder
  wrote past its buffers (AddressSanitizer heap-buffer-overflow). Vector
  lengths, `group >= 1` and `delta` in `0:3` are now validated in C++, and a
  `NA` in a random-effects variable gives a clear R-side error.

* `brs_prep()` Mode 3 rows (only `left`/`right` known) had `y = NA`, so
  `model.frame()` dropped them and those censored observations never reached
  the fit. They are now kept.

## Likelihood: no probability floor, exact tails

* The compiled likelihood floored every censored probability at `1e-15`
  before taking the log, and chose the CDF tail by the position of the
  interval on (0, 1) rather than by the fitted distribution. An interval far
  above a small fitted mean (or below a large one) was computed as
  `F(right) - F(left)` with both terms equal to 1 to machine precision, so it
  hit the floor: the observation contributed the constant `log(1e-15)` with
  zero gradient, and the optimiser maximised a trimmed likelihood that ignored
  outliers. Now the tail is chosen by the mean `a / (a + b)`, `pbeta()` is
  evaluated in plain scale on the small side, there is no floor, and below
  `1e-240` (where R's `bratio` loses accuracy) an endpoint Laplace
  approximation of the tail integral takes over. The endpoint clamp to
  `[1e-5, 1 - 1e-5]` is unchanged. An R mirror of the same rules,
  `.brs_obs_loglik()`, is used by `brs_cv()` for the log-score.

  User-visible consequences: `logLik()`, `AIC()` and `brs_cv()` change on data
  with observations far in a tail, and the precision estimate can drop a lot.
  In the audit example (200 observations, 4 outliers) the estimated precision
  went from 328 to 32; the old value was an artefact of the trimmed
  likelihood. An interval-censored observation with `left == right`
  (probability zero) now contributes `-1e6` instead of `log(1e-15)`.

## Compiled backend: Armadillo, chain-rule derivatives, stable standard errors

* The mixed-model backend (`brsmm()`) is now written in RcppArmadillo, like the
  rest of the package; RcppEigen is no longer a dependency. The port was
  checked against the Eigen code to rounding before any change of method.
* Gradients and Hessians use the chain rule on the linear predictors
  (per-observation central differences, cost independent of the number of
  coefficients). The new default `hessian_method = "cpp"` is about 16x faster
  than `numDeriv` and agrees with it to 1e-8 in the standard errors
  (`hessian_method = "numDeriv"` remains); `brs()` fits are 2.0-2.3x faster.
* `brs()` accepts `start` and `control`. Bootstrap and jackknife refits
  warm-start from the parent estimate and use the compiled Hessian: 1.4-1.8x
  faster for `R = 100`, `n = 250`, with intervals unchanged to 6e-5 standard
  errors.
* `brs_marginaleffects()` draws from the Cholesky factor of the variance
  matrix instead of its eigenvectors, so a negligible change of the variance
  matrix no longer changes the simulated standard errors (an eigenvector sign
  flip moved them by about 2% with `n_sim = 400`).
* `brsmm()` passes the gradient of the chosen approximation (Laplace, AGHQ or
  QMC; chain rule and implicit-function theorem at the modes) to `optim()`
  and computes the Hessian from it (`hessian_method = "cpp"`, default). Standard errors of random-slope models
  are now finite and reproducible (they were `NaN` or changed by 20-40%
  between two practically identical optima). Fits are 1.5-6.6x faster.
* The inner search for the random-effect modes is a Levenberg–Marquardt
  Newton method with warm starts. It no longer returns `b = 0` when the
  curvature there is indefinite, and the silent eigenvalue floor of 1e-8
  (which added up to `+9.2` per direction to the Laplace value) is gone; a
  group without a positive-definite mode is penalised and reported in
  `fit$diagnostics$inner`.
* AGHQ and QMC scale the nodes by the symmetric root `C^(-1/2)` of the
  curvature, so their values no longer depend on the eigenvector sign and order
  conventions of the LAPACK in use (random-effect dimension >= 2). For two or
  more random effects QMC is therefore a different estimator than before; at
  1024 points it underestimates the log-likelihood (mean error -0.05 over 30
  two-effect data sets, -0.03 with the old scaling), so `int_method = "aghq"`
  is recommended up to three random effects. Warm starts of the inner modes
  are reset at every `brsmm()` call, so a fit depends only on its data and
  start.
* Structural errors in the compiled functions (wrong parameter length, `NA` or
  non-finite data, group codes beyond the number of rows) stop with a clear
  message; a `NaN` parameter gives the likelihood penalty, `+-Inf` the bounds.
  An `NA` in a prepared `delta` column now stops in R and in C++; before, it
  reached the compiled code as `NaN` and was cast to an integer (undefined
  behaviour: rejected with a misleading message on x86-64, silently read as an
  exact observation on arm64).

## Documentation

* Help pages aligned with Lopes (2023) and checked against the code: the
  complete likelihood with the four censoring types (the dissertation's
  table swaps `delta = 1` and `2`; the package follows its equation), the
  three parameterisations with their links, mean and variance (`repar = 2`
  dispersion is `1 / (1 + a + b)`, not a coefficient of variation), the
  interval directions, the border handling that replaces the dissertation's
  edge transformation, the residual types, estimation (chain-rule gradient
  and Hessian) and a "Fit diagnostics" section explaining each warning.
* `ncuts` is documented as `K`, the maximum score (scale `0..K`, `K + 1`
  categories); a scale that starts at 1 (Likert 1-5) is shifted to 0-4.
* New help pages for `summary()`, `confint()` and `anova()` (Wald, LRT and
  the chi-bar-square mixture); `brsmm()` documents the integration methods,
  the gradient and the boundary diagnostics.
* New examples, runnable without `\donttest`: a synthetic NRS-11 study
  (times x groups), the four `brs_prep()` modes and a Likert item, one
  simulation per parameterisation, random intercept and slope with
  `anova()`, bootstrap with a factor, cross-validation, marginal effects.
* Vignettes rewritten where they disagreed with the code: likelihood,
  directions, parameterisations, diagnostics, residuals (brs-intro); a Monte
  Carlo study in the dissertation's design (brs-advanced-workflows); ICC,
  SD/Corr intervals and the chi-bar-square test (brs-mm); marginal effects
  and score probabilities (brs-analyst-tools). The intro vignette uses a
  percentile bootstrap instead of BCa, which lowers the build time.
* README: removed the claims of an analytical gradient, of new
  methodology (the parameterisation is Bayer, 2011) and of exactly normal
  quantile residuals; `ncuts` wording fixed; references added. The package
  page attributes the M1/M2/M3 comparison to the dissertation's study.
* `brs_prep()` accepts a bound or score column that is entirely `NA`
  (logical in R) as numeric `NA` instead of stopping.

# betaregscale 2.7.4

Resubmission addressing CRAN feedback on vignette build time (Uwe Ligges,
2026-08-23): *"Please reduce the vignette build timings ... Otherwise we cannot
afford checking the vignette regularly on CRAN."* No change to the statistical
methods or to the user-facing API.

## Vignette build time

* The vignettes now use smaller toy data sets and fewer resampling iterations.
  Total knit time for the four vignettes drops from 77.5s to 14.5s on the
  development machine (5.3x), and `checking re-building of vignette outputs`
  inside `R CMD check --as-cran` drops from 89s to 26s (3.4x). That step took
  376s on the CRAN incoming check, so the expected saving there is around four
  and a half minutes. The source tarball also shrank from 1.35 MB to 1.09 MB.

  The dominant cost was a single chunk in `brs-intro.Rmd`, 45.1s of the 77.5s
  total: a bootstrap with `ci_type = "bca"` on 1000 observations. BCa obtains its
  acceleration constant from a leave-one-out jackknife, so that call performed
  1000 model fits for the jackknife on top of the 100 bootstrap replicates. The
  sample size for that vignette is now 250 and the replicate count 30.

  Other reductions: sample sizes in `brs-intro.Rmd` (1000 to 200/250),
  `brs-advanced-workflows.Rmd` (260 to 150, and the mixed-effects example from
  1200 to 250 observations) and `brs-mm.Rmd` (5 groups of 200 to 12 groups of
  20, which is also a more natural design for illustrating random effects);
  bootstrap replicates (80/100/120 to 30); marginal-effect simulation draws
  (120/160 to 60); and cross-validation repeats (5 to 2). All vignettes still
  knit without warnings.

* `brs-advanced-workflows.Rmd` computed its bootstrap twice (once for the table,
  once again with `ci_type = "bca"` for the forest plot) and its average marginal
  effects twice. Both are now computed once, with `keep_draws = TRUE`, and the
  plots reuse them. Dropping the redundant BCa call also removes its leave-one-out
  jackknife. BCa remains demonstrated in `brs-intro.Rmd`.

## Documentation

* `?brs_bootstrap` gains a section on the cost of `ci_type = "bca"`. The
  leave-one-out jackknife behind the acceleration constant requires `R + n` model
  fits rather than `R`, so the run time is governed by the sample size rather
  than by the number of replicates. This was not documented, and it is easy to
  hit unexpectedly on a large sample.

---

# betaregscale 2.7.3

This release makes no change to the user-facing API. It improves the numerical
conditioning of the interval-censored likelihood, removes an unused C++ backend,
and collects packaging and documentation cleanups for CRAN submission.

## Numerical accuracy

* The interval probability `P(lo < Y < hi)` that underlies every interval-censored
  observation is no longer always computed from lower-tail beta CDF values. When
  both endpoints lie in the upper tail (`lo + hi > 1`, common when the fitted mean
  is close to 1) the difference is now taken between upper-tail (survival)
  probabilities, so both terms stay small and the subtraction no longer suffers
  catastrophic cancellation. The two forms are identical in exact arithmetic; the
  new one is strictly better conditioned in floating point.

## Bug fixes

* Removed two stale help pages, `man/brsmm_loglik_eigen.Rd` and
  `man/brsmm_group_modes_eigen.Rd`, that documented `brsmm_loglik_eigen()` and
  `brsmm_group_modes_eigen()`. Those objects do not exist: the Eigen entry points
  are registered as the internal `.brsmm_loglik_eigen` and
  `.brsmm_group_modes_eigen`. `R CMD check` reported both as code/documentation
  mismatches.
* Removed a Dropbox conflict copy of `.Rbuildignore` that had been committed by
  mistake and was being shipped in the source tarball, where `R CMD check`
  flagged it as a hidden file with a non-portable name.

## Packaging

* `DESCRIPTION` no longer sets `LazyData: true`. The package ships no `data/`
  directory, so `R CMD build` was already reporting
  "Omitted 'LazyData' from DESCRIPTION".
* `betareg` was removed from `Suggests`. It is not used by any function, test or
  vignette; the package is only mentioned in prose when describing the output
  style of `summary()`.
* `src/Makevars` and `src/Makevars.win` no longer request the OpenMP compiler and
  linker flags (`$(SHLIB_OPENMP_CXXFLAGS)`). No translation unit in `src/`
  contains an OpenMP directive, so the flags added portability risk without any
  parallelism.
* `TODO.md`, a development-only file, is now listed in `.Rbuildignore`.

## Documentation

* `NEWS.md` records under 2.7.0 the extension of `int_method = "aghq"` and
  `int_method = "qmc"` to multivariate random effects, which had been implemented
  but never announced. The 2.6.8 heading, which had been concatenated onto the
  end of the 2.6.9 entry, is now a separate section.
* `README.md`: the S3 interface table marked `autoplot()` as unavailable for
  `brsmm` objects although `autoplot.brsmm()` is registered, exported and
  documented; it is now marked as available, and the missing `plot()` and
  extractor (`logLik()`, `AIC()`, `BIC()`, `nobs()`, `fitted()`, `formula()`,
  `model.matrix()`) rows were added.
* `README.md`: the score-probability example described its output as
  "500 patients" while the accompanying simulation creates 1000.
* `README.md`: the installation section claimed the package was "currently under
  review for CRAN". It has been on CRAN since 2.6.9 (published 2026-02-25), so
  `install.packages("betaregscale")` is again presented as the primary route,
  with the GitHub install offered as the development version.
* `README.md`: the package summary and the mixed-effects section described the
  random-effects likelihood as Laplace-only, omitting the AGHQ and QMC methods.
* `README.md`: fixed the interval-censoring notation, which rendered as a
  semicolon-separated list rather than a closed interval.
* The package help page (`?betaregscale`) showed its "Useful links" section
  twice and listed the maintainer a third time below the author list. The cause
  was a block in `R/autoplot.R` that used the `"_PACKAGE"` sentinel purely to
  emit `@rawNamespace` directives, which made roxygen2 treat it as a second
  package-level documentation block; it now uses `@noRd` and contributes only
  the NAMESPACE directives (`NAMESPACE` is unchanged). A hand-written
  `@seealso`/`@author` pair in `R/betaregscale-package.R` that duplicated what
  roxygen2 already derives from `Authors@R`, `URL` and `BugReports` was removed.
* `DESCRIPTION` now lists the GitHub repository in `URL` alongside the pkgdown
  site.
* The URI of the cited master's dissertation is no longer wrapped in `\url{}` in
  the `\references` section of the 22 help pages that carry it. The UFPR
  institutional repository hosting it is down (502 across the whole server, with
  an official maintenance notice), so the link resolved to an error. The
  identifier is persistent and the reference unchanged; only the hyperlink markup
  was dropped, and it will be restored once the repository is back.

## Internal

* Removed the unused Armadillo mixed-effects backend from `src/loglik.cpp`
  (`betaregscale_loglik_mixed_laplace_cpp()`, `betaregscale_group_modes_cpp()`
  and the `build_group_index()` / `group_Q()` / `golden_max_group()` /
  `laplace_group()` helpers, 267 lines). It was reachable from no R code, test or
  vignette: `brsmm()` uses the Eigen backend exclusively. `RcppExports` were
  regenerated and the two corresponding help pages removed.

---

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
  bounds checking.
* The internal `.brsmm_loglik_eigen` entry point is exported with a leading dot to keep it out
  of the public namespace.

---

# betaregscale 2.7.0

## New features

* `brsmm()` now supports the `"aghq"` and `"qmc"` integration methods with
  **multivariate** random effects. Previously both were restricted to a single
  random-effects dimension and `brsmm()` raised an error for
  `random = ~ 1 + x | group` unless `int_method = "laplace"`. The Eigen backend now
  builds a Cartesian adaptive Gauss-Hermite grid and prime-based multidimensional
  Halton sequences of the required dimension, so all three integration methods are
  available for any random-effects structure.
  * For `int_method = "aghq"` the grid has `n_points^q_re` nodes; `brsmm()` stops with
    an informative message if that exceeds 500,000, suggesting a smaller `n_points`
    or `int_method = "qmc"`.
  * `int_method = "qmc"` supports up to 50 random-effects dimensions.
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
* Minor mathematical formatting and typographical fixes (e.g., en-dashes for page ranges) in `README.md` references.

---

# betaregscale 2.6.8

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
