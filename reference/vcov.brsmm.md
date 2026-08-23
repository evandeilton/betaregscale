# Variance-covariance matrix for brsmm coefficients

Variance-covariance matrix for brsmm coefficients

## Usage

``` r
# S3 method for class 'brsmm'
vcov(object, model = c("full", "mean", "precision", "random"), ...)
```

## Arguments

- object:

  A fitted `"brsmm"` object.

- model:

  Character: `"full"`, `"mean"`, `"precision"`, or `"random"`.

- ...:

  Currently ignored.

## Value

Numeric matrix.

## See also

[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md),
[`coef.brsmm`](https://evandeilton.github.io/betaregscale/reference/coef.brsmm.md),
[`confint.brsmm`](https://evandeilton.github.io/betaregscale/reference/confint.brsmm.md)

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
vcov(fit, model = "mean")
#>             (Intercept)         x1
#> (Intercept)   0.6395224 -0.3459754
#> x1           -0.3459754  0.2383237
# }
```
