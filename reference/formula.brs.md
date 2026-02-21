# Extract model formula

Extract model formula

## Usage

``` r
# S3 method for class 'brs'
formula(x, ...)
```

## Arguments

- x:

  A fitted `"betaregscale"` object.

- ...:

  Ignored.

## Value

The formula used to fit the model.

## See also

[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md),
[`model.matrix.brs`](https://evandeilton.github.io/betaregscale/reference/model.matrix.brs.md),
[`coef.brs`](https://evandeilton.github.io/betaregscale/reference/coef.brs.md)

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
formula(fit)
#> y ~ x1
#> <environment: 0x56181a2906a0>
# }
```
