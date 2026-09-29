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
#> -1.9840 -0.4876 -0.1678  0.6964  2.1289 
#> 
#> Coefficients (mean model with logit link):
#>             Estimate Std. Error z value Pr(>|z|)
#> (Intercept)   0.4212     0.7770   0.542    0.588
#> x1           -0.3374     0.4753  -0.710    0.478
#> 
#> Phi coefficients (precision model with logit link):
#>             Estimate Std. Error z value Pr(>|z|)  
#> (Intercept)  -0.5805     0.3229  -1.798   0.0722 .
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Random effects (SD and Corr; 95% Wald CI on the log / atanh scale; no tests, see anova()):
#>                Estimate  Lower  Upper
#> SD (Intercept)   0.5338 0.1259 2.2624
#> ---
#> Mixed beta interval model (Laplace)
#> Observations: 20  | Groups: 4 
#> Log-likelihood: -92.1831 on 4 Df | AIC: 192.3663 | BIC: 196.3492 
#> Pseudo R-squared: 0.0029 
#> Number of iterations: 39 (BFGS) 
#> Censoring: 18 interval | 1 left | 1 right 
#> 
# }
```
