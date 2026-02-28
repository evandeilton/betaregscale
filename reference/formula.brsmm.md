# Extract model formula

Extract model formula

## Usage

``` r
# S3 method for class 'brsmm'
formula(x, ...)
```

## Arguments

- x:

  A fitted `"brsmm"` object.

- ...:

  Ignored.

## Value

The formula used to fit the model.

## See also

[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md),
[`model.matrix.brsmm`](https://evandeilton.github.io/betaregscale/reference/model.matrix.brsmm.md)

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
formula(fit)
#> y ~ x1 | 1
#> <environment: 0x55d7c53feb38>
# }
```
