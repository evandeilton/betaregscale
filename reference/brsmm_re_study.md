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
#> Summary by term:
#>         term sd_model mean_mode sd_mode shrinkage_ratio shapiro_p
#>  (Intercept)   0.5338    0.0011  0.4459          0.6976    0.8241
#> 
#> Estimated covariance matrix D:
#>       [,1]
#> [1,] 0.285
#> 
#> Estimated correlation matrix:
#>      [,1]
#> [1,]    1
rs$summary
#>          term  sd_model  mean_mode   sd_mode shrinkage_ratio shapiro_p
#> 1 (Intercept) 0.5338184 0.00114907 0.4458568       0.6975955 0.8240609
# }
```
