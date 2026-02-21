# Extract fitted values

Extract fitted values

## Usage

``` r
# S3 method for class 'brs'
fitted(object, type = c("mu", "phi"), ...)
```

## Arguments

- object:

  A fitted `"betaregscale"` object.

- type:

  Character: `"mu"` (default) or `"phi"`.

- ...:

  Currently ignored.

## Value

Numeric vector of fitted values.

## See also

[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md),
[`residuals.brs`](https://evandeilton.github.io/betaregscale/reference/residuals.brs.md),
[`predict.brs`](https://evandeilton.github.io/betaregscale/reference/predict.brs.md)

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
head(fitted(fit))
#> [1] 0.5087209 0.4538020 0.5087209 0.4538020 0.5087209 0.4538020
head(fitted(fit, type = "phi"))
#> [1] 0.4030146 0.4030146 0.4030146 0.4030146 0.4030146 0.4030146
# }
```
