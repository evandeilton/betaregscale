## R CMD check results

### Test environments

* Local: Windows 11, R 4.5.2
* GitHub Actions: ubuntu-latest (R release), macOS-latest (R release), windows-latest (R release)

### R CMD check results

0 errors | 0 warnings | 0 notes

### Downstream dependencies

This is a new CRAN submission. There are no downstream dependencies.

### Notes

* This version (2.0.1) adds the `bs_prepare()` data preprocessing
  function and fixes a row-indexing edge case.
* The package includes compiled C++ code via Rcpp and RcppArmadillo.
* All tests pass on all platforms.
