# Extract coefficients from a brsmm fit

Extract coefficients from a brsmm fit

## Usage

``` r
# S3 method for class 'brsmm'
coef(object, model = c("full", "mean", "precision", "random"), ...)
```

## Arguments

- object:

  A fitted `"brsmm"` object.

- model:

  Character: `"full"` (default), `"mean"`, `"precision"`, or `"random"`.

- ...:

  Currently ignored.

## Value

Named numeric vector.

## See also

[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md),
[`vcov.brsmm`](https://evandeilton.github.io/betaregscale/reference/vcov.brsmm.md),
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
coef(fit)
#>                    (Intercept)                             x1 
#>                      0.4213505                     -0.3374768 
#>              (phi)_(Intercept) (re_chol_logsd)_(Intercept)|id 
#>                     -0.5804084                     -0.6277120 
coef(fit, model = "mean")
#> (Intercept)          x1 
#>   0.4213505  -0.3374768 
coef(fit, model = "random")
#> (re_chol_logsd)_(Intercept)|id 
#>                      -0.627712 
# }
```
