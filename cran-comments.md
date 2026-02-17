## R CMD check results

### Test environments

* Local: Ubuntu 24.04, R 4.5.2
* GitHub Actions: ubuntu-latest (R release), macOS-latest (R release), windows-latest (R release)

### R CMD check results

0 errors | 0 warnings | 1 note

### Downstream dependencies

This is a new CRAN submission. There are no downstream dependencies.

### Notes

* This version (2.6.3) revises and expands package vignettes and README with stronger mathematical exposition and analyst-focused workflows.
* Vignettes were fully rendered (`devtools::build_vignettes()`) and documentation site rebuilt (`pkgdown::build_site()`).
* Examples and outputs were standardized with `knitr::kable(..., digits = 4)` for clearer reporting.
* Bibliographic entries in vignettes/README were updated with DOI and dual online validation links.
* The package includes compiled C++ code via Rcpp and RcppArmadillo.
* The only `R CMD check` note is the standard "New submission" incoming
  feasibility note.
