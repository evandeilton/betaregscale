# Fit a variable-dispersion beta interval regression model

Estimates the parameters of a beta regression model with
observation-specific dispersion governed by a second linear predictor.
Both submodels are estimated jointly via maximum likelihood, using the
complete likelihood with mixed censoring.

## Usage

``` r
brs_fit_var(
  formula,
  data,
  link = NULL,
  link_phi = NULL,
  hessian_method = c("numDeriv", "optim"),
  ncuts = NULL,
  lim = NULL,
  repar = 2L,
  method = c("BFGS", "L-BFGS-B")
)
```

## Arguments

- formula:

  A [`Formula`](https://rdrr.io/pkg/Formula/man/Formula.html)-style
  formula with two parts: `y ~ x1 + x2 | z1 + z2`.

- data:

  Data frame.

- link:

  Link for the first parameter (the mean under `repar = 1, 2`; the shape
  \\p\\ under `repar = 0`). `NULL` (default) selects the link implied by
  `repar`: `"logit"` for `repar = 1, 2`, `"log"` for `repar = 0`. See
  the 'Reparameterizations and links' section of
  [`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md)
  for the admissible values.

- link_phi:

  Link for the second parameter. `NULL` (default) selects `"logit"` for
  `repar = 2` (dispersion on \\(0, 1)\\) and `"log"` for `repar = 0, 1`
  (positive shape/precision).

- hessian_method:

  Character: `"numDeriv"` or `"optim"`.

- ncuts:

  Number of scale categories. `NULL` (default) uses the value stored by
  [`brs_prep`](https://evandeilton.github.io/betaregscale/reference/brs_prep.md)
  in `attr(data, "ncuts")`, or 100 when `data` was not prepared. A value
  that differs from the stored one is ignored with a warning (the
  endpoints were built with the stored value).

- lim:

  Uncertainty half-width. `NULL` (default) uses `attr(data, "lim")` from
  [`brs_prep`](https://evandeilton.github.io/betaregscale/reference/brs_prep.md),
  or 0.5; same rule as `ncuts`.

- repar:

  Reparameterization scheme (default 2); see
  [`brs_repar`](https://evandeilton.github.io/betaregscale/reference/brs_repar.md).

- method:

  Optimization method (default `"BFGS"`).

## Value

An object of class `"brs"`.

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
fit <- brs_fit_var(y ~ x1 | x2, data = prep)
print(fit)
#> 
#> Call:
#> brs_fit_var(formula = y ~ x1 | x2, data = prep)
#> 
#> Coefficients (mean model with logit link):
#> (Intercept)          x1 
#>      0.2732     -0.2310 
#> 
#> Phi coefficients (precision model with logit link):
#> (Intercept)          x2 
#>     -0.3789     -0.0288 
#> 
# }
```
