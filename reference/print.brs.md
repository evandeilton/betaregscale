# Print a fitted model (brief betareg style)

Print a fitted model (brief betareg style)

## Usage

``` r
# S3 method for class 'brs'
print(x, digits = max(3, getOption("digits") - 3), ...)
```

## Arguments

- x:

  A fitted `"betaregscale"` object.

- digits:

  Number of significant digits.

- ...:

  Included for consistency with generic methods. Currently passed to
  internal methods where applicable.

## Value

Invisibly returns the input object `x`. The function is called for its
side effect of printing a formatted summary of the fitted model to the
console, including the model call, mean coefficients (with link
function), and precision coefficients (with link function).

## See also

[`summary.brs`](https://evandeilton.github.io/betaregscale/reference/summary.brs.md),
[`print.summary.brs`](https://evandeilton.github.io/betaregscale/reference/print.summary.brs.md),
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md)

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
print(fit)
#> 
#> Call:
#> brs(formula = y ~ x1, data = prep)
#> 
#> Coefficients (mean model with logit link):
#> (Intercept)          x1 
#>      0.2551     -0.2202 
#> 
#> Phi coefficients (precision model with logit link):
#>   (phi) 
#> -0.3929 
#> 
# }
```
