## Resubmission (2.6.7)

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

### Test environments

* Local: Ubuntu 24.04, R 4.5.2
* GitHub Actions: ubuntu-latest (R release), macOS-latest (R release), windows-latest (R release)

### R CMD check results

0 errors | 0 warnings | 2 notes

### Downstream dependencies

There are no downstream dependencies.

### Notes

* Version 2.6.6 addresses the two CRAN notes above (DESCRIPTION wording and URL).
* The package includes compiled C++ code via Rcpp and RcppArmadillo.
* R CMD check notes observed locally:
  - "New submission" (incoming feasibility);
  - "unable to verify current time" in `checking for future file timestamps`.
