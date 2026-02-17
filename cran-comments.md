## R CMD check results

### Test environments

* Local: Ubuntu 24.04, R 4.5.2
* GitHub Actions: ubuntu-latest (R release), macOS-latest (R release), windows-latest (R release)

### R CMD check results

0 errors | 0 warnings | 2 notes

### Downstream dependencies

This is a new CRAN submission. There are no downstream dependencies.

### Notes

* This version (2.6.5) improves analyst-facing uncertainty and interpretation workflows for `brs` models:
  - `brs_bootstrap()` now supports BCa intervals (`ci_type = "bca"`), endpoint Monte Carlo diagnostics (`mcse_lower`, `mcse_upper`), and side-by-side Wald interval columns (`wald_lower`, `wald_upper`).
  - `autoplot.brs_bootstrap()` now supports visual comparison of bootstrap vs Wald intervals in `type = "ci_forest"`.
  - `brs_marginaleffects()` robustness was improved (central differences, adaptive step scaling, and reusable simulation draws), and `autoplot.brs_marginaleffects()` was added for forest/magnitude/distribution summaries.
  - `brs_cens()` was refined with richer summaries and optional, domain-agnostic interpretation messages (`inform = TRUE`) suitable for any scale-censored context.
* README and vignettes were updated to document the new bootstrap and marginal-effects workflows, including visual diagnostics and reproducible examples.
* Bibliographic entries in package documentation continue to use DOI-based validated links (`https://doi.org/...`) for CRAN URL checks.
* The package includes compiled C++ code via Rcpp and RcppArmadillo.
* `R CMD check` notes are:
  - standard "New submission" incoming feasibility note;
  - "unable to verify current time" in `checking for future file timestamps`,
    observed in this local environment.
