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

## Details

`icc` is the intraclass correlation of \\\mathrm{logit}(Y)\\ implied by
the fitted model: for two observations of the same group with the
covariates of observation \\i\\, \$\$\mathrm{ICC}\_i =
\frac{\mathrm{Var}\_b\[\psi(a_i) - \psi(b_i)\]}
{\mathrm{Var}\_b\[\psi(a_i) - \psi(b_i)\] + E_b\[\psi_1(a_i) +
\psi_1(b_i)\]},\$\$ where \\a_i(b), b_i(b)\\ are the beta shapes with
random part \\b \sim N(0, x\_{r,i}^\top D x\_{r,i})\\, and \\\psi\\,
\\\psi_1\\ are the digamma and trigamma functions
(\\E\[\mathrm{logit}\\Y\] = \psi(a) - \psi(b)\\,
\\\mathrm{Var}\[\mathrm{logit}\\Y\] = \psi_1(a) + \psi_1(b)\\). The
expectations over \\b\\ use 40-point Gauss-Hermite quadrature and the
reported value is the mean of \\\mathrm{ICC}\_i\\ over the observations.
The level-1 variance is that of the beta, so the value depends on the
precision; with the logit link (only) and a large precision it
approaches \\\sigma_b^2 / (\sigma_b^2 + \psi_1(a) + \psi_1(b))\\. It
replaces the logistic-latent formula \\\sigma_b^2 / (\sigma_b^2 +
\pi^2/3)\\, which does not describe a beta response.

The moments of \\\mathrm{logit}(Y)\\ over \\b\\ can be infinite: with a
probit link when \\\sigma_b^2 \ge 1/2\\, with a cloglog link for every
\\\sigma_b \> 0\\, and numerically with any link when \\\sigma_b\\ is
very large. The value would then be set by the clamp of the mean
(\\10^{-5}\\), so `icc` is `NA`, with a warning, whenever
\\E_b\[\psi_1(a) + \psi_1(b)\]\\ changes by more than 10% when that
clamp is tightened to \\10^{-8}\\.

## References

Lopes, J. E. (2023). *Modelos de regressao beta para dados de escala*.
Master's dissertation, Universidade Federal do Parana, Curitiba. URI:
https://hdl.handle.net/1884/86624.

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
#>   (Intercept)                 0.5339
#> 
#> ICC (logit(Y) scale, beta level-1 variance): 0.1807
#> 
#> Summary by term (SD_model = model SD; shrinkage = Var(modes)/Var(model)):
#>         term sd_model mean_mode sd_mode shrinkage_ratio shapiro_p
#>  (Intercept)   0.5339    0.0012   0.446          0.6976     0.824
rs$summary
#>          term  sd_model   mean_mode   sd_mode shrinkage_ratio shapiro_p
#> 1 (Intercept) 0.5339184 0.001150958 0.4459553       0.6976425 0.8240296
# }
```
