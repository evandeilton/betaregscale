# Fit a fixed-dispersion beta interval regression model

Estimates the parameters of a beta regression model with a single
(scalar) dispersion parameter using maximum likelihood. The
log-likelihood and its gradient are evaluated by the compiled C++
backend supporting the complete likelihood with mixed censoring types.

## Usage

``` r
brs_fit_fixed(
  formula,
  data,
  link = NULL,
  link_phi = NULL,
  ncuts = NULL,
  lim = NULL,
  hessian_method = c("cpp", "numDeriv", "optim"),
  repar = 2L,
  method = c("BFGS", "L-BFGS-B"),
  interval = NULL,
  start = NULL,
  control = list()
)
```

## Arguments

- formula:

  Two-sided formula `y ~ x1 + x2 + ...`.

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

- ncuts:

  Integer \\K\\, the maximum score: the scale is \\0, 1, \ldots, K\\
  (\\K + 1\\ categories). `NULL` (default) uses the value stored by
  [`brs_prep`](https://evandeilton.github.io/betaregscale/reference/brs_prep.md)
  in `attr(data, "ncuts")`, or 100 when `data` was not prepared. A value
  that differs from the stored one is ignored with a warning (the
  endpoints were built with the stored value).

- lim:

  Half-width of the score cell in \\(0, 0.5\]\\ (`interval = "mid"`
  only). `NULL` (default) uses `attr(data, "lim")` from
  [`brs_prep`](https://evandeilton.github.io/betaregscale/reference/brs_prep.md),
  or 0.5; same rule as `ncuts`. Values below 0.5 warn (partial
  coarsening).

- hessian_method:

  Character: `"cpp"` (default), `"numDeriv"` or `"optim"`. `"cpp"` uses
  the compiled chain-rule Hessian (per-observation second derivatives in
  the linear predictors); `"numDeriv"` differentiates the log-likelihood
  with [`hessian`](https://rdrr.io/pkg/numDeriv/man/hessian.html);
  `"optim"` keeps the optimizer's own approximation.

- repar:

  Reparameterization scheme (default 2); see
  [`brs_repar`](https://evandeilton.github.io/betaregscale/reference/brs_repar.md).

- method:

  Optimization method: `"BFGS"` (default) or `"L-BFGS-B"`.

- interval:

  Direction of the uncertainty interval, `"mid"`, `"right"` or `"left"`
  (see
  [`brs_check`](https://evandeilton.github.io/betaregscale/reference/brs_check.md)).
  `NULL` (default) uses `attr(data, "interval")` from
  [`brs_prep`](https://evandeilton.github.io/betaregscale/reference/brs_prep.md),
  or `"mid"`; same rule as `ncuts`.

- start:

  Optional numeric vector of starting values (mean coefficients, then
  dispersion coefficients). `NULL` (default) uses
  [`compute_start()`](https://evandeilton.github.io/betaregscale/reference/compute_start.md).
  Refits (bootstrap, jackknife) pass the parent estimate here.

- control:

  Control list for [`optim`](https://rdrr.io/r/stats/optim.html); its
  entries are merged into the default `list(maxit = 5000L)`.

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
fit <- brs_fit_fixed(y ~ x1 + x2, data = prep)
print(fit)
#> 
#> Call:
#> brs_fit_fixed(formula = y ~ x1 + x2, data = prep)
#> 
#> Coefficients (mean model with logit link):
#> (Intercept)          x1          x2 
#>      0.0664     -0.2268      0.3884 
#> 
#> Phi coefficients (precision model with logit link):
#>   (phi) 
#> -0.4091 
#> 
# }
```
