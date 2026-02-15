## R CMD check results

### Test environments

* Local: Windows 11, R 4.5.2
* GitHub Actions: ubuntu-latest (R release), macOS-latest (R release), windows-latest (R release)

### R CMD check results

0 errors | 0 warnings | 1 note

### Downstream dependencies

This is a new CRAN submission. There are no downstream dependencies.

### Notes

* This version (2.4.0) consolidates simulation into `brs_sim()`,
  removes `brs_sim_var()` from the public API, and keeps
  `brs_loglik()`/`brs_loglik_var()` as internal helpers.
* The package includes compiled C++ code via Rcpp and RcppArmadillo.
* The only `R CMD check` note is the standard "New submission" incoming
  feasibility note.
