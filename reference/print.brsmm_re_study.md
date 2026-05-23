# Print a random-effects study

Prints a compact summary of the random-effects study returned by
[`brsmm_re_study`](https://evandeilton.github.io/betaregscale/reference/brsmm_re_study.md),
including per-term standard deviations, shrinkage ratios, Shapiro-Wilk
p-values, and the estimated covariance and correlation matrices.

## Usage

``` r
# S3 method for class 'brsmm_re_study'
print(x, digits = max(3, getOption("digits") - 3), ...)
```

## Arguments

- x:

  A `"brsmm_re_study"` object returned by
  [`brsmm_re_study`](https://evandeilton.github.io/betaregscale/reference/brsmm_re_study.md).

- digits:

  Integer: number of significant digits for rounding (default
  `max(3, getOption("digits") - 3)`).

- ...:

  Currently ignored.

## Value

Invisibly returns `x`. Called for its side-effect of printing the study
to the console.

## See also

[`brsmm_re_study`](https://evandeilton.github.io/betaregscale/reference/brsmm_re_study.md),
[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md),
[`ranef.brsmm`](https://evandeilton.github.io/betaregscale/reference/ranef.brsmm.md)

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
rs <- brsmm_re_study(fit)
print(rs)
#> 
#> Random-effects study
#> Groups: 4 
#> 
#> Random-effects (VarCorr):
#>   Name                      Std.Dev.
#>   re1                         0.5340
#> 
#> ICC (latent logistic scale): 0.0798
#> 
#> Summary by term (SD_model = model SD; shrinkage = Var(modes)/Var(model)):
#>         term sd_model mean_mode sd_mode shrinkage_ratio shapiro_p
#>  (Intercept)    0.534     0.001  0.4461          0.6977    0.8239
# }
```
