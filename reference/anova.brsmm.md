# Likelihood-ratio comparison involving mixed models

[`anova()`](https://rdrr.io/r/stats/anova.html) for `"brsmm"` fits: the
same table as
[`anova.brs`](https://evandeilton.github.io/betaregscale/reference/anova.brs.md),
with the chi-bar-square mixture \\\frac12\chi^2\_{df-1} +
\frac12\chi^2\_{df}\\ for rows that add one random-effect term, whose
variance is on the boundary under \\H_0\\. This is the test to use for a
variance component: the Wald statistic of its log standard deviation is
not meaningful (see
[`summary.brsmm`](https://evandeilton.github.io/betaregscale/reference/summary.brsmm.md)).

## Usage

``` r
# S3 method for class 'brsmm'
anova(object, ..., test = c("Chisq", "none"))
```

## Arguments

- object:

  A fitted `"brsmm"` model.

- ...:

  Further fitted `"brsmm"` and/or `"brs"` models.

- test:

  `"Chisq"` (default) or `"none"`.

## Value

An object of class `"anova"`; see
[`anova.brs`](https://evandeilton.github.io/betaregscale/reference/anova.brs.md).

## References

Self, S. G., and Liang, K.-Y. (1987). Asymptotic properties of maximum
likelihood estimators and likelihood ratio tests under nonstandard
conditions. *Journal of the American Statistical Association*,
**82**(398), 605–610.
[doi:10.1080/01621459.1987.10478472](https://doi.org/10.1080/01621459.1987.10478472)

## See also

[`anova.brs`](https://evandeilton.github.io/betaregscale/reference/anova.brs.md),
[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md),
[`summary.brsmm`](https://evandeilton.github.io/betaregscale/reference/summary.brsmm.md)

## Examples

``` r
set.seed(11)
g <- 20
d <- data.frame(id = factor(rep(1:g, each = 8)), x = runif(8 * g))
shp <- brs_repar(plogis(-0.4 + d$x + rnorm(g, sd = 0.6)[d$id]), phi = 0.25)
d$y <- round(10 * rbeta(nrow(d), shp$shape1, shp$shape2))
m0 <- brs(y ~ x, data = d, ncuts = 10)
m1 <- brsmm(y ~ x, random = ~ 1 | id, data = d, ncuts = 10)
anova(m0, m1)  # Pr(>Chisq) = half the chi2(1) tail
#> Likelihood-ratio comparison of brs/brsmm models
#> Rows M2: one added random effect (variance on the boundary); Pr(>Chisq) from the chi-bar-square mixture 1/2 chi2(Df - 1) + 1/2 chi2(Df).
#> 
#>            Df  logLik    AIC    BIC  Chisq Chi Df Pr(>Chisq)    
#> M1 (brs)    3 -370.88 747.76 756.99                             
#> M2 (brsmm)  4 -363.94 735.88 748.18 13.883      1  9.725e-05 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
