# Akaike information criterion

\\-2\ell + k\\p\\, with \\p\\ the number of estimated parameters.

## Usage

``` r
# S3 method for class 'brs'
AIC(object, ..., k = 2)
```

## Arguments

- object:

  A fitted `"brs"` object.

- ...:

  Further fitted models (`brs`, `brsmm` or others with an
  [`AIC()`](https://rdrr.io/r/stats/AIC.html) method).

- k:

  Penalty per parameter (default 2).

## Value

Scalar AIC value; with several models, a data frame with columns `df`
and `AIC`, one row per model, as
[`stats::AIC()`](https://rdrr.io/r/stats/AIC.html).

## See also

[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md),
[`logLik.brs`](https://evandeilton.github.io/betaregscale/reference/logLik.brs.md),
[`BIC.brs`](https://evandeilton.github.io/betaregscale/reference/BIC.brs.md),
[`brs_gof`](https://evandeilton.github.io/betaregscale/reference/brs_gof.md)

## Examples

``` r
# \donttest{
dat <- data.frame(
  y = c(
    0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
    10, 40, 55, 70, 85, 25, 35, 65, 80, 15
  ),
  x1 = rep(c(1, 2), 10)
)
prep <- brs_prep(dat, ncuts = 100)
#> brs_prep: n = 20 | exact = 0, left = 1, right = 1, interval = 18
fit <- brs(y ~ x1, data = prep)
AIC(fit)
#> [1] 191.3041
fit0 <- brs(y ~ 1, data = prep)
AIC(fit0, fit)
#>      df      AIC
#> fit0  2 189.4692
#> fit   3 191.3041
# }
```
