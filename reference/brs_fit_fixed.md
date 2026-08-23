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
  link = "logit",
  link_phi = "logit",
  ncuts = 100L,
  lim = 0.5,
  hessian_method = c("numDeriv", "optim"),
  repar = 2L,
  method = c("BFGS", "L-BFGS-B")
)
```

## Arguments

- formula:

  Two-sided formula `y ~ x1 + x2 + ...`.

- data:

  Data frame.

- link:

  Mean link function (default `"logit"`).

- link_phi:

  Dispersion link function (default `"logit"`).

- ncuts:

  Number of scale categories (default 100).

- lim:

  Uncertainty half-width (default 0.5).

- hessian_method:

  Character: `"numDeriv"` (default) or `"optim"`. With `"numDeriv"` the
  Hessian is computed after convergence using
  [`hessian`](https://rdrr.io/pkg/numDeriv/man/hessian.html), which is
  typically more accurate than the built-in optim Hessian.

- repar:

  Reparameterization scheme (default 2).

- method:

  Optimization method: `"BFGS"` (default) or `"L-BFGS-B"`.

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
