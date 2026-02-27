# Variance-covariance matrix of estimated coefficients

Variance-covariance matrix of estimated coefficients

## Usage

``` r
# S3 method for class 'brs'
vcov(object, model = c("full", "mean", "precision"), ...)
```

## Arguments

- object:

  A fitted `"betaregscale"` object.

- model:

  Character: which component (`"full"`, `"mean"`, or `"precision"`).

- ...:

  Ignored.

## Value

A square numeric matrix.

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
#> (Intercept)  0.74717128 -0.444368041 -0.010878706
#> x1          -0.44436804  0.292891277  0.009448033
#> (phi)       -0.01087871  0.009448033  0.076320255
vcov(fit, model = "mean")
#>             (Intercept)         x1
#> (Intercept)   0.7471713 -0.4443680
#> x1           -0.4443680  0.2928913
# }
```
