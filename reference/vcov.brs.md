# Variance-covariance matrix of estimated coefficients

Variance-covariance matrix of estimated coefficients

## Usage

``` r
# S3 method for class 'brs'
vcov(object, model = c("full", "mean", "precision"), ...)
```

## Arguments

- object:

  A fitted `"brs"` object.

- model:

  Character: which component (`"full"`, `"mean"`, or `"precision"`).

- ...:

  Ignored.

## Value

A square numeric matrix.

## Details

\\(-H)^{-1}\\ with \\H\\ the Hessian of the log-likelihood at the
estimate. No generalised inverse is used: a singular Hessian gives an
`NA` matrix, and negative or non-finite variances become `NA` (row and
column); both cases warn (see `fit$diagnostics`).

## See also

[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md),
[`coef.brs`](https://evandeilton.github.io/betaregscale/reference/coef.brs.md),
[`confint.brs`](https://evandeilton.github.io/betaregscale/reference/confint.brs.md)

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
vcov(fit)
#>             (Intercept)           x1        (phi)
#> (Intercept)  0.74717309 -0.444369079 -0.010878513
#> x1          -0.44436908  0.292891939  0.009448082
#> (phi)       -0.01087851  0.009448082  0.076320569
vcov(fit, model = "mean")
#>             (Intercept)         x1
#> (Intercept)   0.7471731 -0.4443691
#> x1           -0.4443691  0.2928919
# }
```
