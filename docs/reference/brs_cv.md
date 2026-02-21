# K-fold cross-validation for brs models

Performs repeated k-fold cross-validation for
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md)
models.

## Usage

``` r
brs_cv(formula, data, k = 5L, repeats = 1L, ...)
```

## Arguments

- formula:

  Model formula passed to
  [`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md).

- data:

  Data frame.

- k:

  Number of folds.

- repeats:

  Number of repeated k-fold runs.

- ...:

  Additional arguments forwarded to
  [`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md)
  (e.g., `repar`, `link`, `method`).

## Value

A data frame with one row per fold and columns: `repeat`, `fold`,
`n_train`, `n_test`, `log_score`, `rmse_yt`, `mae_yt`, `converged`, and
`error`. The object has class `"brs_cv"`.

## Details

The `log_score` is the mean log predictive contribution under the
complete likelihood contribution implied by each observation's censoring
type (`delta`).

## References

Lopes, J. E. (2023). *Modelos de regressao beta para dados de escala*.
Master's dissertation, Universidade Federal do Parana, Curitiba. URI:
<https://hdl.handle.net/1884/86624>.

Hawker, G. A., Mian, S., Kendzerska, T., and French, M. (2011). Measures
of adult pain: Visual Analog Scale for Pain (VAS Pain), Numeric Rating
Scale for Pain (NRS Pain), McGill Pain Questionnaire (MPQ), Short-Form
McGill Pain Questionnaire (SF-MPQ), Chronic Pain Grade Scale (CPGS),
Short Form-36 Bodily Pain Scale (SF-36 BPS), and Measure of Intermittent
and Constant Osteoarthritis Pain (ICOAP). Arthritis Care and Research,
63(S11), S240-S252.
[doi:10.1002/acr.20543](https://doi.org/10.1002/acr.20543)

Hjermstad, M. J., Fayers, P. M., Haugen, D. F., et al. (2011). Studies
comparing Numerical Rating Scales, Verbal Rating Scales, and Visual
Analogue Scales for assessment of pain intensity in adults: a systematic
literature review. Journal of Pain and Symptom Management, 41(6),
1073-1093.
[doi:10.1016/j.jpainsymman.2010.08.016](https://doi.org/10.1016/j.jpainsymman.2010.08.016)

## Examples

``` r
# \donttest{
dat <- data.frame(
  y = c(
    0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
    10, 40, 55, 70, 85, 25, 35, 65, 80, 15
  ),
  x1 = rep(c(1, 2), 10),
  x2 = rep(c(0, 0, 1, 1), 5)
)
prep <- brs_prep(dat, ncuts = 100)
#> brs_prep: n = 20 | exact = 0, left = 1, right = 1, interval = 18
cv <- brs_cv(y ~ x1, data = prep, k = 3, repeats = 1)
cv
#>   repeat fold n_train n_test log_score rmse_yt mae_yt converged
#> 1      1    1      13      7        NA      NA     NA     FALSE
#> 2      1    2      13      7        NA      NA     NA     FALSE
#> 3      1    3      14      6        NA      NA     NA     FALSE
#>                                   error
#> 1 missing value where TRUE/FALSE needed
#> 2 missing value where TRUE/FALSE needed
#> 3 missing value where TRUE/FALSE needed
# }
```
