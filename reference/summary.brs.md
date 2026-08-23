# Summarize a fitted model (betareg style)

Summarize a fitted model (betareg style)

## Usage

``` r
# S3 method for class 'brs'
summary(object, ...)
```

## Arguments

- object:

  A fitted `"betaregscale"` object.

- ...:

  Ignored.

## Value

A list of class `"summary.betaregscale"`.

## See also

[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md),
[`print.summary.brs`](https://evandeilton.github.io/betaregscale/reference/print.summary.brs.md),
[`brs_est`](https://evandeilton.github.io/betaregscale/reference/brs_est.md),
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
s <- summary(fit)
s$coefficients$mean
#>              Estimate Std. Error    z value  Pr(>|z|)
#> (Intercept)  0.255100  0.8643918  0.2951208 0.7679016
#> x1          -0.220206  0.5411949 -0.4068885 0.6840899
# }
```
