# Compare fitted brs models in a single table

Builds a comparison table for one or more fitted `"brs"` objects,
summarizing fit statistics and (optionally) censoring composition.

## Usage

``` r
brs_table(
  ...,
  models = NULL,
  include_censoring = TRUE,
  sort_by = c("none", "AIC", "BIC", "logLik"),
  decreasing = FALSE,
  digits = 4L
)
```

## Arguments

- ...:

  Fitted `"brs"` objects passed individually.

- models:

  Optional list of fitted `"brs"` objects. Use either `...` or `models`,
  not both.

- include_censoring:

  Logical; include censoring counts/proportions. Default is `TRUE`.

- sort_by:

  Character; optional sort criterion: `"none"` (default), `"AIC"`,
  `"BIC"`, or `"logLik"`.

- decreasing:

  Logical; sort direction when `sort_by != "none"`.

- digits:

  Integer number of digits used for numeric rounding.

## Value

A data frame with one row per model.

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
m1 <- brs(y ~ 1, data = prep)
#> Warning: the standard deviation is zero
m2 <- brs(y ~ x1, data = prep)
brs_table(null = m1, x1 = m2, sort_by = "AIC")
#>   model nobs npar   logLik      AIC      BIC pseudo_r2 exact left right
#> 1  null   20    2 -92.7346 189.4692 191.4606        NA     0    1     1
#> 2    x1   20    3 -92.6521 191.3041 194.2913    0.0029     0    1     1
#>   interval prop_exact prop_left prop_right prop_interval
#> 1       18          0      0.05       0.05           0.9
#> 2       18          0      0.05       0.05           0.9
# }
```
