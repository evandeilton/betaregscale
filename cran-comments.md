## R CMD check results

### Test environments

* Local: Ubuntu 24.04, R 4.5.2
* GitHub Actions: ubuntu-latest (R release), macOS-latest (R release), windows-latest (R release)

### R CMD check results

0 errors | 0 warnings | 1 note

### Downstream dependencies

This is a new CRAN submission. There are no downstream dependencies.

### Notes

* This version (2.6.4) extends `brsmm()` to multivariate random effects in the mean predictor (random intercept + random slope), adds generic `anova()` model comparison methods (`anova.brs`, `anova.brsmm`), and introduces dedicated random-effects studies (`brsmm_re_study`) with new visual diagnostics.
* Vignettes were updated and rebuilt with an explicit evolutionary model-selection workflow (`brs` -> `brsmm RI` -> `brsmm RI+RS`) based on likelihood-ratio comparisons.
* Bibliographic entries in package documentation were revised to keep only DOI-based validated links (`https://doi.org/...`) for CRAN URL checks.
* The package includes compiled C++ code via Rcpp and RcppArmadillo.
* The only `R CMD check` note is the standard "New submission" incoming
  feasibility note.
