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
  (e.g., `repar`, `link`, `interval`, `method`).

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
https://hdl.handle.net/1884/86624.

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
# Synthetic NRS-11 scores: 4 groups x 3 times. Simulated, not real data.
set.seed(3)
nrs <- expand.grid(id = 1:40, time = c("6h", "12h", "24h"))
nrs$group <- factor(paste0("g", (nrs$id - 1) %% 4 + 1))
eta <- -1.3 + c(0, 0.75, 0.3)[nrs$time] + c(0, -0.1, 0.05, 0.1)[nrs$group]
shp <- brs_repar(mu = plogis(eta), phi = 0.3, repar = 2)
nrs$y <- round(10 * rbeta(nrow(nrs), shp$shape1, shp$shape2))

# 3-fold CV of two nested models, time only (m1) and time + group (m2);
# log_score is the mean held-out log-likelihood contribution (higher is better)
set.seed(5)
cv1 <- brs_cv(y ~ time, data = nrs, k = 3, ncuts = 10)
set.seed(5)
cv2 <- brs_cv(y ~ time + group, data = nrs, k = 3, ncuts = 10)
c(m1 = mean(cv1$log_score), m2 = mean(cv2$log_score))
#>        m1        m2 
#> -2.157497 -2.190544 
cv2
#>   repeat fold n_train n_test log_score   rmse_yt    mae_yt converged error
#> 1      1    1      80     40 -2.073336 0.2592196 0.2187621      TRUE  <NA>
#> 2      1    2      80     40 -2.256562 0.3123485 0.2750955      TRUE  <NA>
#> 3      1    3      80     40 -2.241734 0.2697334 0.2188253      TRUE  <NA>
```
