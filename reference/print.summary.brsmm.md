# Print summary for brsmm models

Print summary for brsmm models

## Usage

``` r
# S3 method for class 'summary.brsmm'
print(x, digits = max(3, getOption("digits") - 3), ...)
```

## Arguments

- x:

  A `"summary.brsmm"` object.

- digits:

  Number of digits.

- ...:

  Passed to `printCoefmat`.

## Value

Invisibly returns `x`.

## See also

[`summary.brsmm`](https://evandeilton.github.io/betaregscale/reference/summary.brsmm.md),
[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md),
[`print.brsmm`](https://evandeilton.github.io/betaregscale/reference/print.brsmm.md)

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
print(summary(fit))
#> 
#> Call:
#> brsmm(formula = y ~ x1, random = ~1 | id, data = prep)
#> 
#> Randomized Quantile Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -2.7260 -0.4934 -0.1681  0.6929  2.0620 
#> 
#> Coefficients (mean model with logit link):
#>             Estimate Std. Error z value Pr(>|z|)
#> (Intercept)   0.4213     0.8803   0.479    0.632
#> x1           -0.3374     0.5361  -0.629    0.529
#> 
#> Phi coefficients (precision model with logit link):
#>             Estimate Std. Error z value Pr(>|z|)  
#> (Intercept)  -0.5806     0.3281   -1.77   0.0768 .
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Random effects (SD and Corr; 95% Wald CI on the log / atanh scale; no tests, see anova()):
#>                Estimate  Lower  Upper
#> SD (Intercept)   0.5339 0.1206 2.3642
#> ---
#> Mixed beta interval model (Laplace)
#> Observations: 20  | Groups: 4 
#> Log-likelihood: -92.1831 on 4 Df | AIC: 192.3663 | BIC: 196.3492 
#> Pseudo R-squared: 0.0029 
#> Number of iterations: 28 (BFGS) 
#> Censoring: 18 interval | 1 left | 1 right 
#> 
# }
```
