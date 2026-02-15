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
set.seed(42)
n <- 120
dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n), z1 = rnorm(n))

sim <- brs_sim(
  formula = ~ x1 + x2 | z1, data = dat,
  beta = c(0.2, -0.4, 0.2), zeta = c(0.1, -0.2),
  ncuts = 100, repar = 2
)

m1 <- brs(y ~ x1 + x2, data = sim, repar = 2)
m2 <- brs(y ~ x1 + x2 | z1, data = sim, repar = 2)

brs_table(fixed = m1, variable = m2, sort_by = "AIC")
#>      model nobs npar    logLik      AIC      BIC pseudo_r2 exact left right
#> 1 variable  120    5 -508.6481 1027.296 1041.234    0.1142     0    8    12
#> 2    fixed  120    4 -511.1839 1030.368 1041.518    0.1143     0    8    12
#>   interval prop_exact prop_left prop_right prop_interval
#> 1      100          0    0.0667        0.1        0.8333
#> 2      100          0    0.0667        0.1        0.8333
# }
```
