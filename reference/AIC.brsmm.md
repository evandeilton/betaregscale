# AIC for brsmm models

\\-2\ell + k\\p\\, with \\p\\ the number of parameters (fixed effects,
precision, packed random-effect and cutpoint parameters).

## Usage

``` r
# S3 method for class 'brsmm'
AIC(object, ..., k = 2)
```

## Arguments

- object:

  A fitted `"brsmm"` object.

- ...:

  Further fitted models (`brsmm`, `brs` or others with an
  [`AIC()`](https://rdrr.io/r/stats/AIC.html) method).

- k:

  Numeric penalty per parameter (default 2).

## Value

Numeric scalar; with several models, a data frame with columns `df` and
`AIC`, one row per model, as
[`stats::AIC()`](https://rdrr.io/r/stats/AIC.html).

## See also

[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md),
[`logLik.brsmm`](https://evandeilton.github.io/betaregscale/reference/logLik.brsmm.md),
[`BIC.brsmm`](https://evandeilton.github.io/betaregscale/reference/BIC.brsmm.md),
[`brs_gof`](https://evandeilton.github.io/betaregscale/reference/brs_gof.md)

## Examples

``` r
# \donttest{
dat <- data.frame(
  y = c(
    0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
    10, 40, 55, 70, 85, 25, 35, 65, 80, 15
  ),
  x1 = rep(c(1, 2), 10),
  id = factor(rep(1:4, each = 5))
)
prep <- brs_prep(dat, ncuts = 100)
#> brs_prep: n = 20 | exact = 0, left = 1, right = 1, interval = 18
fit <- brsmm(y ~ x1, random = ~ 1 | id, data = prep)
AIC(fit)
#> [1] 192.3663
AIC(brs(y ~ x1, data = prep), fit)
#>                          df      AIC
#> brs(y ~ x1, data = prep)  3 191.3041
#> fit                       4 192.3663
# }
```
