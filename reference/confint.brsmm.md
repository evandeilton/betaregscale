# Wald confidence intervals for brsmm models

Wald confidence intervals for brsmm models

## Usage

``` r
# S3 method for class 'brsmm'
confint(
  object,
  parm,
  level = 0.95,
  model = c("full", "mean", "precision", "random"),
  ...
)
```

## Arguments

- object:

  A fitted `"brsmm"` object.

- parm:

  Character or integer: which parameters.

- level:

  Confidence level (default 0.95).

- model:

  Character: `"full"`, `"mean"`, `"precision"`, or `"random"`.

- ...:

  Currently ignored.

## Value

Matrix with columns for lower and upper confidence bounds.

## See also

[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md),
[`coef.brsmm`](https://evandeilton.github.io/betaregscale/reference/coef.brsmm.md),
[`vcov.brsmm`](https://evandeilton.github.io/betaregscale/reference/vcov.brsmm.md)

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
confint(fit, model = "mean")
#>                 2.5 %    97.5 %
#> (Intercept) -1.174693 2.0168622
#> x1          -1.313782 0.6391828
# }
```
