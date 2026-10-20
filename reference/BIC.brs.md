# Bayesian information criterion

Bayesian information criterion

## Usage

``` r
# S3 method for class 'brs'
BIC(object, ..., n = c("obs", "groups"))
```

## Arguments

- object:

  A fitted `"brs"` object.

- ...:

  Further fitted models (`brs`, `brsmm` or others with a
  [`BIC()`](https://rdrr.io/r/stats/AIC.html) method).

- n:

  `"obs"` (default): the number of observations. `"groups"` applies to
  `brsmm` fits only (see
  [`BIC.brsmm`](https://evandeilton.github.io/betaregscale/reference/BIC.brsmm.md))
  and is an error here.

## Value

Scalar BIC value; with several models, a data frame with columns `df`
and `BIC`, one row per model, as
[`stats::BIC()`](https://rdrr.io/r/stats/AIC.html).

## See also

[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md),
[`logLik.brs`](https://evandeilton.github.io/betaregscale/reference/logLik.brs.md),
[`AIC.brs`](https://evandeilton.github.io/betaregscale/reference/AIC.brs.md),
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
BIC(fit)
#> [1] 194.2913
# }
```
