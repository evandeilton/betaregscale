# Residuals from a brsmm model

Residuals from a brsmm model

## Usage

``` r
# S3 method for class 'brsmm'
residuals(
  object,
  type = c("response", "pearson", "deviance", "rqr", "weighted", "sweighted"),
  ...
)
```

## Arguments

- object:

  A fitted `"brsmm"` object.

- type:

  Character: `"response"` (default), `"pearson"`, `"deviance"`, `"rqr"`,
  `"weighted"`, or `"sweighted"`.

- ...:

  Currently ignored.

## Value

Numeric vector.

## See also

[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md),
[`fitted.brsmm`](https://evandeilton.github.io/betaregscale/reference/fitted.brsmm.md),
[`plot.brsmm`](https://evandeilton.github.io/betaregscale/reference/plot.brsmm.md)

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
head(residuals(fit))
#> [1] -0.3855877 -0.2593508 -0.1855977  0.1906492  0.3644023  0.3286546
head(residuals(fit, type = "pearson"))
#> [1] -1.3225313 -0.9367197 -0.6365835  0.6885841  1.2498673  1.1086975
# }
```
