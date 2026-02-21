# BIC for brsmm models

BIC for brsmm models

## Usage

``` r
# S3 method for class 'brsmm'
BIC(object, ...)
```

## Arguments

- object:

  A fitted `"brsmm"` object.

- ...:

  Currently ignored.

## Value

Numeric scalar.

## See also

[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md),
[`logLik.brsmm`](https://evandeilton.github.io/betaregscale/reference/logLik.brsmm.md),
[`AIC.brsmm`](https://evandeilton.github.io/betaregscale/reference/AIC.brsmm.md),
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
BIC(fit)
#> [1] 196.3492
# }
```
