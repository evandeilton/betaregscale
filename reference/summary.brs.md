# Summarize a fitted beta interval model

Wald tables for the mean and precision (or dispersion) coefficients,
information criteria, a pseudo \\R^2\\, the censoring counts and a
summary of the randomized quantile residuals.

## Usage

``` r
# S3 method for class 'brs'
summary(object, ...)
```

## Arguments

- object:

  A fitted `"brs"` object.

- ...:

  Currently ignored.

## Value

A list of class `"summary.brs"` with `coefficients` (tables `mean` and
`precision` with columns `Estimate`, `Std. Error`, `z value`,
`Pr(>|z|)`), `residuals` (RQR), `loglik`, `AIC`, `BIC`, `df`, `nobs`,
`pseudo.r2`, `censoring` (counts by type), `link`, `link_phi`, `repar`,
`convergence` and `iterations`.

## Details

For each coefficient \\\hat\theta_j\\, the standard error is \\SE_j =
\sqrt{\[(-H)^{-1}\]\_{jj}}\\ from
[`vcov.brs`](https://evandeilton.github.io/betaregscale/reference/vcov.brs.md),
the Wald statistic is \\z_j = \hat\theta_j / SE_j\\ and the two-sided
p-value is \\2\Phi(-\|z_j\|)\\ (Lopes, 2023, "Inferencia"). The test is
on the link scale (\\H_0: \theta_j = 0\\). When \\-H\\ is singular, or
its inverse has negative variances, the affected standard errors,
statistics and p-values are `NA` (see 'Fit diagnostics' in
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md)).

\\\mathrm{AIC} = -2\ell + 2k\\ and \\\mathrm{BIC} = -2\ell + k\log n\\,
with \\\ell\\ the maximised log-likelihood, \\k\\ the number of
coefficients and \\n\\ the number of observations. The pseudo \\R^2\\ is
the squared correlation between the fitted linear predictor of the mean
and \\g_1(y_i)\\ at the cell centres `yt` (Ferrari and Cribari-Neto,
2004); under `repar = 0` both sides are on the logit scale. It uses the
cell centres, so it is rough when most observations are censored (the
print says so).

The randomized quantile residuals
([`residuals.brs`](https://evandeilton.github.io/betaregscale/reference/residuals.brs.md),
`type = "rqr"`) are drawn without changing the caller's RNG state.

## References

Lopes, J. E. (2023). *Modelos de regressao beta para dados de escala*.
Master's dissertation, Universidade Federal do Parana, Curitiba. URI:
https://hdl.handle.net/1884/86624.

Ferrari, S. L. P., and Cribari-Neto, F. (2004). Beta regression for
modelling rates and proportions. *Journal of Applied Statistics*,
**31**(7), 799–815.
[doi:10.1080/0266476042000214501](https://doi.org/10.1080/0266476042000214501)

## See also

[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md),
[`confint.brs`](https://evandeilton.github.io/betaregscale/reference/confint.brs.md),
[`anova.brs`](https://evandeilton.github.io/betaregscale/reference/anova.brs.md),
[`brs_gof`](https://evandeilton.github.io/betaregscale/reference/brs_gof.md)

## Examples

``` r
set.seed(2023)
d <- data.frame(time = factor(rep(c("6h", "12h", "24h"), each = 60),
                              levels = c("6h", "12h", "24h")))
shp <- brs_repar(mu = plogis(-1.3 + c(0, 0.75, 0.3)[d$time]), phi = 0.3)
d$y <- round(10 * rbeta(nrow(d), shp$shape1, shp$shape2))
fit <- brs(y ~ time, data = d, ncuts = 10)
s <- summary(fit)
s
#> 
#> Call:
#> brs(formula = y ~ time, data = d, ncuts = 10)
#> 
#> Quantile residuals:
#>     Min      1Q  Median      3Q     Max 
#> -2.3714 -0.7045 -0.0964  0.6605  2.6798 
#> 
#> Coefficients (mean model with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  -1.3079     0.1628  -8.036 9.31e-16 ***
#> time12h       0.5245     0.2156   2.433   0.0150 *  
#> time24h       0.4992     0.2151   2.321   0.0203 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Phi coefficients (precision model with logit link):
#>       Estimate Std. Error z value Pr(>|z|)    
#> (phi)  -0.8393     0.1103   -7.61 2.73e-14 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> ---
#> Log-likelihood: -377.6945 on 4 Df | AIC: 763.3891 | BIC: 776.1609 
#> Pseudo R-squared: 0.0537  (midpoint approx.; interpret with caution for heavily censored data) 
#> Number of iterations: 26 (BFGS) 
#> Censoring: 141 interval | 37 left | 2 right 
#> 
s$coefficients$mean
#>               Estimate Std. Error   z value     Pr(>|z|)
#> (Intercept) -1.3079316  0.1627663 -8.035643 9.308939e-16
#> time12h      0.5245223  0.2155579  2.433324 1.496090e-02
#> time24h      0.4992069  0.2150869  2.320955 2.028928e-02
c(AIC = s$AIC, BIC = s$BIC, pseudo_R2 = s$pseudo.r2)
#>          AIC          BIC    pseudo_R2 
#> 763.38908156 776.16090897   0.05368292 
```
