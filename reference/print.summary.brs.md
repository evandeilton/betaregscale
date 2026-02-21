# Print a model summary (betareg style)

Print a model summary (betareg style)

## Usage

``` r
# S3 method for class 'summary.brs'
print(x, digits = max(3, getOption("digits") - 3), ...)
```

## Arguments

- x:

  A `"summary.betaregscale"` object.

- digits:

  Number of digits.

- ...:

  Passed to `printCoefmat`.

## Value

Invisibly returns the input object `x`. The function is called for its
side effect of printing a comprehensive summary to the console,
including the model call, quantile residuals, coefficient tables for
mean and precision submodels with significance stars, goodness-of-fit
statistics (log-likelihood, pseudo R-squared), optimization details, and
censoring information.

## See also

[`summary.brs`](https://evandeilton.github.io/betaregscale/reference/summary.brs.md),
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md),
[`print.brs`](https://evandeilton.github.io/betaregscale/reference/print.brs.md)

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
print(summary(fit))
#> 
#> Call:
#> brs(formula = y ~ x1, data = prep)
#> 
#> Quantile residuals:
#>     Min      1Q  Median      3Q     Max 
#> -2.6418 -0.4830  0.0678  0.5471  2.1625 
#> 
#> Coefficients (mean model with logit link):
#>             Estimate Std. Error z value Pr(>|z|)
#> (Intercept)   0.2551     0.8644   0.295    0.768
#> x1           -0.2202     0.5412  -0.407    0.684
#> 
#> Phi coefficients (precision model with logit link):
#>       Estimate Std. Error z value Pr(>|z|)
#> (phi)  -0.3929     0.2763  -1.422    0.155
#> ---
#> Log-likelihood: -92.6521 on 3 Df | AIC: 191.3041 | BIC: 194.2913 
#> Pseudo R-squared: 0.0029 
#> Number of iterations: 17 (BFGS) 
#> Censoring: 18 interval | 1 left | 1 right 
#> 
# }
```
