# Random-effects study for brsmm models

Provides a compact numeric study of random effects, including: estimated
covariance matrix, correlation matrix, per-term standard deviations,
empirical mean/SD of posterior modes, shrinkage ratio, and a normality
check by Shapiro-Wilk (when applicable).

## Usage

``` r
brsmm_re_study(object, ...)
```

## Arguments

- object:

  A fitted `"brsmm"` object.

- ...:

  Currently ignored.

## Value

A list with class `"brsmm_re_study"`.

## References

Lopes, J. E. (2023). *Modelos de regressao beta para dados de escala*.
Master's dissertation, Universidade Federal do Parana, Curitiba. URI:
<https://hdl.handle.net/1884/86624>.

Ferrari, S. L. P., and Cribari-Neto, F. (2004). Beta regression for
modelling rates and proportions. *Journal of Applied Statistics*,
**31**(7), 799–815.
[doi:10.1080/0266476042000214501](https://doi.org/10.1080/0266476042000214501)

## See also

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
rs$summary
#>          term  sd_model   mean_mode   sd_mode shrinkage_ratio shapiro_p
#> 1 (Intercept) 0.5340369 0.001048123 0.4460754       0.6977084 0.8239009
# }
```
