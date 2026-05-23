# Extract model coefficients

Extract model coefficients

## Usage

``` r
# S3 method for class 'brs'
coef(object, model = c("full", "mean", "precision"), ...)
```

## Arguments

- object:

  A fitted `"betaregscale"` object.

- model:

  Character: which component to return. `"full"` (default) returns all
  parameters, `"mean"` returns only the mean-model coefficients,
  `"precision"` returns only the precision coefficients.

- ...:

  Ignored.

## Value

Named numeric vector of estimated parameters.

## See also

[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md),
[`brs_est`](https://evandeilton.github.io/betaregscale/reference/brs_est.md),
[`vcov.brs`](https://evandeilton.github.io/betaregscale/reference/vcov.brs.md)

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
coef(fit)
#> (Intercept)          x1       (phi) 
#>   0.2551000  -0.2202060  -0.3929144 
coef(fit, model = "mean")
#> (Intercept)          x1 
#>    0.255100   -0.220206 
coef(fit, model = "precision")
#>      (phi) 
#> -0.3929144 
# }
```
