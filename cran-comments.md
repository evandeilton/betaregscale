## Submission (2.7.3)

This submission supersedes 2.7.1. There is no change to the user-facing API.
It improves the numerical conditioning of the interval-censored likelihood,
removes an unused C++ backend, and collects packaging cleanups and corrections
to the package documentation found while preparing the release.

### Numerical accuracy

1. **Interval probability computed on the better-conditioned tail**
   (`src/brs_common.h`). Every interval-censored observation contributes
   `log P(lo < Y < hi)`, which was always evaluated as
   `pbeta(hi, lower) - pbeta(lo, lower)`. When both endpoints sit in the upper
   tail the two lower-tail values are both close to 1 and the subtraction loses
   most of its significant digits. The difference is now taken between
   upper-tail (survival) probabilities whenever `lo + hi > 1`, so both terms
   remain small. The two expressions are identical in exact arithmetic.

   Concretely, for `mu = 0.95`, `phi = 1000` (shapes `a = 950`, `b = 50`) and
   the interval `[0.99, 0.995]`, the old form returns exactly `0` in double
   precision while the new one returns `7.97e-20`. Across a grid of shapes and
   intervals typical of saturating rating scales, relative discrepancies reach
   `3e-3` for probabilities well above the internal `1e-15` floor, i.e. for
   values that enter the log-likelihood directly.

### Dead code

2. **Unused Armadillo mixed-effects backend removed** (`src/loglik.cpp`, 267
   lines): `betaregscale_loglik_mixed_laplace_cpp()`,
   `betaregscale_group_modes_cpp()` and the `build_group_index()` / `group_Q()` /
   `golden_max_group()` / `laplace_group()` helpers were reachable from no R
   code, test or vignette. `brsmm()` uses the Eigen backend exclusively.
   `RcppExports` were regenerated and the two orphaned help pages removed.

### Corrections found by R CMD check

3. **Code/documentation mismatch (was a WARNING).** `man/brsmm_loglik_eigen.Rd`
   and `man/brsmm_group_modes_eigen.Rd` documented `brsmm_loglik_eigen()` and
   `brsmm_group_modes_eigen()`, which do not exist: the Eigen entry points are
   registered as the internal `.brsmm_loglik_eigen` and
   `.brsmm_group_modes_eigen`. Both files were stale and have been removed.

4. **Non-portable file name (was a WARNING) and hidden file (was a NOTE).** A
   Dropbox conflict copy of `.Rbuildignore` had been committed by mistake and
   was being shipped in the source tarball. It has been removed.

### Packaging

5. **`LazyData` removed from DESCRIPTION.** The package ships no `data/`
   directory, so `R CMD build` was reporting
   "Omitted 'LazyData' from DESCRIPTION" on every build.

6. **`betareg` removed from `Suggests`.** It is not used by any function, test
   or vignette; the package is referenced only in prose when describing the
   output style of `summary()`.

7. **OpenMP flags removed from `src/Makevars` and `src/Makevars.win`.** No
   translation unit in `src/` contains an OpenMP directive, so
   `$(SHLIB_OPENMP_CXXFLAGS)` was requested in `PKG_CXXFLAGS` and `PKG_LIBS`
   without any parallel region to justify it.

8. **`TODO.md` added to `.Rbuildignore`**, as a development-only file.

9. **`URL`** now lists the GitHub repository alongside the pkgdown site.

### Documentation

10. **Duplicated sections on the package help page.** A block in `R/autoplot.R`
   used the `"_PACKAGE"` sentinel although it exists only to emit
   `@rawNamespace` directives for the `ggplot2::autoplot` S3 methods. roxygen2
   treated it as a second package-level documentation block, so
   `?betaregscale` rendered its "Useful links" section twice. The block now
   uses `@noRd`; `NAMESPACE` is byte-identical before and after the change.
   A hand-written `@seealso`/`@author` pair in `R/betaregscale-package.R` that
   repeated what roxygen2 already derives from `Authors@R`, `URL` and
   `BugReports` was removed.

11. **`NEWS.md`.** The extension of `int_method = "aghq"` and
   `int_method = "qmc"` to multivariate random effects was implemented but
   never recorded; it is now documented under 2.7.0. The "betaregscale 2.6.8"
   heading had been concatenated onto the final line of the 2.6.9 entry and so
   never rendered as a section heading.

12. **`README.md`.** The S3 interface table marked `autoplot()` as unavailable
   for `brsmm` objects although `autoplot.brsmm()` is registered, exported and
   documented; the missing `plot()` and extractor rows were added. The
   score-probability example described its output as "500 patients" while the
   accompanying simulation creates 1000. The installation section presented
   `install.packages("betaregscale")` as the stable route while the package is
   still under review. The package summary described the random-effects
   likelihood as Laplace-only, omitting the AGHQ and QMC methods.

---

## Previous submission (2.7.0)

This is a feature release adding diagnostic and plotting enhancements to both the
`brs` (fixed-effects) and `brsmm` (mixed-effects) model classes.

### Changes in this version:

1. **`theme` argument in all autoplot functions:** `autoplot.brs()` and `autoplot.brsmm()`
   now accept a `theme` argument (default `ggplot2::theme_minimal()`) that is propagated
   to every internal plotting helper, replacing the previously hardcoded theme. Accepts
   both theme objects and theme functions.

2. **Label overrides (`title`, `xlab`, `ylab`) and `type = "all"`:** Both autoplot
   dispatchers gain `title`, `xlab`, and `ylab` for post-hoc label overrides via
   `ggplot2::labs()`, a `type = "all"` option that arranges all panels in a single
   `gridExtra::grid.arrange()` grid, and `ncol` to control grid layout. The `...`
   argument forwards named arguments to `ggplot2::theme()` on top of the base theme.

3. **New `type = "shrinkage"` in `autoplot.brsmm()`:** Scatter of Laplace posterior
   modes vs. naïve per-group logit-mean deviations, with identity line and loess smoother.

4. **Improved `type = "ranef_caterpillar"`:** Error bars now represent ±1.96 × marginal
   Model SD (from the `D` covariance matrix), with an informative subtitle.

5. **New `plot.brsmm()` panels:** `which = 7` (Q-Q normal of RE modes) and `which = 8`
   (dotchart caterpillar), both using base R graphics.

6. **ICC in `brsmm_re_study()`:** Returns `$icc` (intraclass correlation on the latent
   logistic scale). `print.brsmm_re_study()` shows a VarCorr-style table and the ICC.

7. **Cosmetic coefficient printing:** `print.summary.*` and `print.*` methods strip
   internal prefixes (`(phi)_`, `(re_chol_logsd)_`, `(re_chol)_`) for cleaner output.

8. **Deprecated `size` aesthetic fix:** Calibration plots now use `linewidth` instead
   of `size` for `geom_line()`, eliminating the ggplot2 ≥ 3.4.0 deprecation warning.

9. **New test file:** `test-re-and-autoplot-improvements.R` with 15 new focused tests.

---

## Previous submission (2.6.9)

This is a minor resubmission to address documentation completeness and language consistency across the package.

### Changes in this version:

1. **Documentation completeness:** Added missing `\value`, `\seealso`, and runnable `\donttest{}` examples to several recently standardized S3 methods (`print.summary`, `residuals`, `summary`, `vcov`, `ranef`) to guarantee full CRAN compliance.
2. **Vignette and README formatting:** Translated remaining Portuguese text to English in the mixed-effects vignette (`brs-mm.Rmd`), corrected generic method dispatch documentation (such as using `ranef()` instead of `ranef.brsmm()`), and fixed typographical issues (e.g., en-dashes in references).

---

## Previous resubmission (2.6.8)

### Changes in this version:

1. **brsmm standardization:** Added missing extractors (`formula()`, `model.matrix()`, `confint()`) and augmented conditional outputs for `residuals()` and `predict()` to achieve explicit method parity with the base `brs` object implementations. Function parameter consistency guarantees full compliance with CRAN policies (no modifications to `.GlobalEnv` or base `par()` properties were made).

2. **Documentation audit — examples:**
   - Added `\donttest{}` examples to all ~30 previously undocumented exported functions (S3 methods for `brs` and `brsmm`: `coef()`, `vcov()`, `logLik()`, `AIC()`, `BIC()`, `nobs()`, `formula()`, `model.matrix()`, `fitted()`, `residuals()`, `confint()`, `predict()`, `print()`, `summary()`, `ranef()`, `anova()`, `plot()`, `autoplot()`).
   - All examples use deterministic toy datasets — no `set.seed()` in any example.
   - No `\dontrun{}` anywhere in the package.

3. **Documentation audit — references:**
   - Added DOI `10.1080/0266476042000214501` (Ferrari & Cribari-Neto, 2004) to every occurrence of that reference.

4. **Documentation audit — structure:**
   - Added `@seealso` cross-links to all S3 method blocks.
   - Completed `@param`, `@return`, `@method`, and `@examples` for `print.brsmm_re_study()`, `ranef()` generic, `brs_coef()`, `brs_hessian()`, and all extracting methods like `vcov()`, `predict()`, etc.
   - `autoplot.*` examples use `ggplot2::autoplot()` (explicit namespace) for reliability.

---

## Previous resubmission (2.6.7)
This is a resubmission following CRAN feedback (Konstanze Lauseker, 20 Feb 2026):

### Changes in response to CRAN review:

1. **Missing \value tags:** Added \value documentation to `man/print.brs.Rd` and `man/print.summary.brs.Rd` explaining that both functions invisibly return the input object and describing the printed output structure.

2. **\dontrun{} usage:** Replaced `\dontrun{}` with `\donttest{}` in `man/brs_gof.Rd` and created a complete executable example including data simulation and model fitting.

3. **.GlobalEnv modification:** Removed all code that modified `.GlobalEnv` from `R/bootstrap.R` (lines 94-110). The `seed` parameter was removed from `brs_bootstrap()`. Users should now call `set.seed()` before the function if reproducibility is needed.

4. **set.seed() in functions:** Removed `set.seed()` calls from four exported functions:
   - `R/bootstrap.R`: `brs_bootstrap()` (removed `seed` parameter)
   - `R/marginaleffects.R`: `brs_marginaleffects()` (removed `seed` parameter)
   - `R/brsmm.R`: `brsmm()` (removed `seed` parameter)
   - `R/cv.R`: `brs_cv()` (removed `seed` parameter)

   All functions now expect users to call `set.seed()` externally for reproducibility. Examples have been updated to demonstrate this pattern.

### Previous resubmission notes (18 Feb 2026):

1. **Function names in DESCRIPTION:** All function names in the Description field now use parentheses (e.g. `logLik()`, `coef()`, `confint()`, `vcov()`). The same convention has been applied in the package-level documentation (`man/betaregscale-package.Rd`).

2. **Invalid URL (404):** The URL was updated from https://evandeilton.github.io/betaregscale/ to the package repository https://github.com/evandeilton/betaregscale in both DESCRIPTION and the package Rd.

---

## R CMD check results

0 errors | 0 warnings | 2 notes

Neither note is a property of the package as submitted:

1. `checking CRAN incoming feasibility ... NOTE` — *Found the following
   (possibly) invalid URLs: `https://hdl.handle.net/1884/86624`, Status: 500*.
   This is the persistent handle of the master's dissertation cited in the
   `\references` section of 22 help pages. The handle resolver itself
   (`https://hdl.handle.net/`) answers 200; the institutional repository it
   redirects to, `https://acervodigital.ufpr.br/`, is currently returning 502
   for **every** URL including its own home page, i.e. the whole server is
   down. The handle is a persistent identifier and the reference is correct; we
   will re-verify immediately before submitting.

2. `checking HTML version of manual ... NOTE` — *Skipping checking HTML
   validation: no command 'tidy' found*. HTML Tidy is not installed on the
   local check machine. This note does not arise on the CRAN check farm.

For reference, the previous state of the package (2.7.1) produced
**2 warnings and 3 notes**; the two warnings (code/documentation mismatch and
non-portable file name) and the two packaging notes (hidden file, non-standard
top-level file) are resolved by this release.

`checking installed package size` reports 22.6Mb (`libs` 20.2Mb) as INFO. The
Ubuntu R build injects `-g` into `CXXFLAGS` (`R CMD config CXXFLAGS`), so
`betaregscale.so` ships unstripped debug information; the same object stripped
with `strip --strip-debug` is 325 KB. `src/Makevars` and `src/Makevars.win`
contain only `CXX_STD`, `-DARMA_64BIT_WORD` and the LAPACK/BLAS/FLIBS link
line. `checking compilation flags in Makevars`, `checking compilation flags
used`, `checking compiled code` and `checking C++ specification` all pass.
`CXX_STD = CXX17` is retained deliberately: the package declares
`Depends: R (>= 4.1.0)`, and R only defaults to C++17 from 4.3.0 onwards.

### Timings

* `checking whether package can be installed`: 33s
* `checking examples` (with `--run-donttest`): 13s
* `checking tests`: 14s (all testthat tests pass)
* `checking re-building of vignette outputs`: 90s

## Test environments

* Local: Ubuntu (Linux 7.0.0-30-generic), R 4.6.1 (2026-06-24),
  x86_64-pc-linux-gnu, `R CMD check --as-cran`
* GitHub Actions: ubuntu-latest (R devel, release, oldrel-1),
  macOS-latest (R release), windows-latest (R release)

## Downstream dependencies

There are no downstream dependencies.

## Notes

* The package contains compiled C++ code via Rcpp, RcppArmadillo and RcppEigen.
* Expected on a first submission: the "New submission" note from the incoming
  feasibility check.
