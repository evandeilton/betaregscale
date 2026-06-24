# Print a fitted brsmm model

Print a fitted brsmm model

## Usage

``` r
# S3 method for class 'brsmm'
print(x, digits = max(3, getOption("digits") - 3), ...)
```

## Arguments

- x:

  A fitted `"brsmm"` object.

- digits:

  Number of digits.

- ...:

  Included for consistency with generic methods.

## Value

Invisibly returns `x`.

## See also

[`summary.brsmm`](https://evandeilton.github.io/betaregscale/reference/summary.brsmm.md),
[`print.summary.brsmm`](https://evandeilton.github.io/betaregscale/reference/print.summary.brsmm.md),
[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md)

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
print(fit)
#> 
#> Call:
#> brsmm(formula = y ~ x1, random = ~1 | id, data = prep)
#> 
#> Coefficients (mean model with logit link):
#> (Intercept)          x1 
#>      0.4214     -0.3375 
#> 
#> Phi coefficients (precision model with logit link):
#> (Intercept) 
#>     -0.5804 
#> 
#> Random-effects parameters:
#> logSD.(Intercept)|id 
#>              -0.6277 
#> 
#> Random SD: 0.5338 
#> ---
#> Mixed beta interval model (Laplace)
#> Observations: 20  | Groups: 4 
#> Log-likelihood: -92.1831 
#> Convergence code: 0 
# }
```
