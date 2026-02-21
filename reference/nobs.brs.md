# Number of observations

Number of observations

## Usage

``` r
# S3 method for class 'brs'
nobs(object, ...)
```

## Arguments

- object:

  A fitted `"betaregscale"` object.

- ...:

  Ignored.

## Value

Integer: number of observations.

## See also

[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md),
[`fitted.brs`](https://evandeilton.github.io/betaregscale/reference/fitted.brs.md),
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
nobs(fit)
#> [1] 20
# }
```
