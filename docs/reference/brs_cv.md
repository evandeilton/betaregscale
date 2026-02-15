# K-fold cross-validation for brs models

Performs repeated k-fold cross-validation for
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md)
models.

## Usage

``` r
brs_cv(formula, data, k = 5L, repeats = 1L, seed = NULL, ...)
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

- seed:

  Optional random seed for reproducibility.

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

Hawker, G. A., Mian, S., Kendzerska, T., and French, M. (2011). Measures
of adult pain: Visual Analog Scale for Pain (VAS Pain), Numeric Rating
Scale for Pain (NRS Pain), McGill Pain Questionnaire (MPQ), Short-Form
McGill Pain Questionnaire (SF-MPQ), Chronic Pain Grade Scale (CPGS),
Short Form-36 Bodily Pain Scale (SF-36 BPS), and Measure of Intermittent
and Constant Osteoarthritis Pain (ICOAP). Arthritis Care and Research,
63(S11), S240-S252. doi:10.1002/acr.20543.

Hjermstad, M. J., Fayers, P. M., Haugen, D. F., et al. (2011). Studies
comparing Numerical Rating Scales, Verbal Rating Scales, and Visual
Analogue Scales for assessment of pain intensity in adults: a systematic
literature review. Journal of Pain and Symptom Management, 41(6),
1073-1093. doi:10.1016/j.jpainsymman.2010.08.016.

## Examples

``` r
# \donttest{
set.seed(99)
d <- data.frame(x1 = rnorm(150), x2 = rnorm(150))
sim <- brs_sim(
  formula = ~ x1 + x2, data = d,
  beta = c(0.2, -0.5, 0.3), phi = 0.2, ncuts = 100, repar = 2
)

cv <- brs_cv(y ~ x1 + x2, data = sim, k = 3, repeats = 1, repar = 2)
cv
#>   repeat fold n_train n_test log_score   rmse_yt    mae_yt converged error
#> 1      1    1     100     50 -4.102422 0.3993795 0.3618229      TRUE  <NA>
#> 2      1    2     100     50 -4.436561 0.3532272 0.3135323      TRUE  <NA>
#> 3      1    3     100     50 -4.011085 0.3558378 0.3308460      TRUE  <NA>
# }
```
