# Fitted values from a brsmm model

Fitted values from a brsmm model

## Usage

``` r
# S3 method for class 'brsmm'
fitted(object, type = c("mu", "phi"), ...)
```

## Arguments

- object:

  A fitted `"brsmm"` object.

- type:

  Character: `"mu"` (default) or `"phi"`.

- ...:

  Currently ignored.

## Value

Numeric vector.

## See also

[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md),
[`residuals.brsmm`](https://evandeilton.github.io/betaregscale/reference/residuals.brsmm.md),
[`predict.brsmm`](https://evandeilton.github.io/betaregscale/reference/predict.brsmm.md)

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
head(fitted(fit))
#> [1] 0.3856147 0.3093481 0.3856147 0.3093481 0.3856147 0.5713120
head(fitted(fit, type = "phi"))
#> [1] 0.3588126 0.3588126 0.3588126 0.3588126 0.3588126 0.3588126
# }
```
