# Wald confidence intervals

Computes Wald confidence intervals for model parameters using the normal
approximation.

## Usage

``` r
# S3 method for class 'brs'
confint(
  object,
  parm,
  level = 0.95,
  model = c("full", "mean", "precision"),
  ...
)
```

## Arguments

- object:

  A fitted `"betaregscale"` object.

- parm:

  Character or integer: which parameters. If missing, all parameters are
  returned.

- level:

  Confidence level (default 0.95).

- model:

  Character: `"full"`, `"mean"`, or `"precision"`.

- ...:

  Currently ignored.

## Value

Matrix with columns for lower and upper confidence bounds.

## See also

[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md),
[`coef.brs`](https://evandeilton.github.io/betaregscale/reference/coef.brs.md),
[`vcov.brs`](https://evandeilton.github.io/betaregscale/reference/vcov.brs.md),
[`brs_est`](https://evandeilton.github.io/betaregscale/reference/brs_est.md)

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
confint(fit)
#>                  2.5 %    97.5 %
#> (Intercept) -1.4390767 1.9492766
#> x1          -1.2809285 0.8405165
#> (phi)       -0.9343775 0.1485488
confint(fit, model = "mean")
#>                 2.5 %    97.5 %
#> (Intercept) -1.439077 1.9492766
#> x1          -1.280929 0.8405165
# }
```
