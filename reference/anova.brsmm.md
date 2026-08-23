# Model comparison by analysis of deviance (LR test) for \`brsmm\`

Model comparison by analysis of deviance (LR test) for \`brsmm\`

## Usage

``` r
# S3 method for class 'brsmm'
anova(object, ..., test = c("Chisq", "none"))
```

## Arguments

- object:

  A fitted `"brsmm"` model.

- ...:

  Additional fitted `"brsmm"` and/or `"brs"` models to compare.

- test:

  Character; `"Chisq"` (default) or `"none"`.

## Value

An object of class `"anova"` and `"data.frame"` with model-wise
log-likelihood, information criteria, and (optionally) LR test columns.

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
[`logLik.brsmm`](https://evandeilton.github.io/betaregscale/reference/logLik.brsmm.md),
[`AIC.brsmm`](https://evandeilton.github.io/betaregscale/reference/AIC.brsmm.md),
[`BIC.brsmm`](https://evandeilton.github.io/betaregscale/reference/BIC.brsmm.md)

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
m1 <- brs(y ~ 1, data = prep)
#> Warning: the standard deviation is zero
m2 <- brsmm(y ~ x1, random = ~ 1 | id, data = prep)
anova(m1, m2)
#>            Df  logLik    AIC    BIC  Chisq Chi Df Pr(>Chisq)
#> M1 (brs)    2 -92.735 189.47 191.46                         
#> M2 (brsmm)  4 -92.183 192.37 196.35 1.1029      2     0.5761
# }
```
