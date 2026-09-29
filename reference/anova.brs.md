# Likelihood-ratio comparison of nested beta interval models

Compares fitted `"brs"` and `"brsmm"` models by log-likelihood, AIC, BIC
and likelihood-ratio tests (Lopes, 2023, "Inferencia").

## Usage

``` r
# S3 method for class 'brs'
anova(object, ..., test = c("Chisq", "none"))
```

## Arguments

- object:

  A fitted `"brs"` model.

- ...:

  Further fitted `"brs"` and/or `"brsmm"` models.

- test:

  `"Chisq"` (default) or `"none"`.

## Value

An object of class `"anova"` (a data frame) with columns `Df`, `logLik`,
`AIC`, `BIC` and, for `test = "Chisq"`, `Chisq`, `Chi Df` and
`Pr(>Chisq)`; the attribute `"heading"` explains the p-values.

## Details

The models are sorted by their number of parameters. For consecutive
models the statistic is \\LR = 2(\ell_1 - \ell_0)\\, with `Chi Df` the
difference in the number of parameters, and \\p = P(\chi^2\_{df} \>
LR)\\. The models must be nested; this is not checked. They must also
describe the same response: the same observations, `interval`, `ncuts`
and (under `"mid"`) `lim`, or the call stops, since different
coarsenings are different likelihoods.

When the larger model adds one random-effect term (a `"brs"` model
against a random-intercept `"brsmm"`, or one more correlated random
term), its variance lies on the boundary of the parameter space under
\\H_0\\ and \\LR\\ follows the mixture \\\frac12\chi^2\_{df-1} +
\frac12\chi^2\_{df}\\ (Self and Liang, 1987; Stram and Lee, 1994); for
one variance component alone this is \\\frac12\chi^2_0 +
\frac12\chi^2_1\\, i.e. half the naive p-value. The printed heading
names the rows that use it. When more than one random term is added at
once, the naive \\\chi^2\_{df}\\ p-value is kept and flagged as
conservative.

## References

Lopes, J. E. (2023). *Modelos de regressao beta para dados de escala*.
Master's dissertation, Universidade Federal do Parana, Curitiba. URI:
https://hdl.handle.net/1884/86624.

Self, S. G., and Liang, K.-Y. (1987). Asymptotic properties of maximum
likelihood estimators and likelihood ratio tests under nonstandard
conditions. *Journal of the American Statistical Association*,
**82**(398), 605–610.
[doi:10.1080/01621459.1987.10478472](https://doi.org/10.1080/01621459.1987.10478472)

Stram, D. O., and Lee, J. W. (1994). Variance components testing in the
longitudinal mixed effects model. *Biometrics*, **50**(4), 1171–1177.
[doi:10.2307/2533455](https://doi.org/10.2307/2533455)

## See also

[`anova.brsmm`](https://evandeilton.github.io/betaregscale/reference/anova.brsmm.md),
[`summary.brs`](https://evandeilton.github.io/betaregscale/reference/summary.brs.md),
[`logLik.brs`](https://evandeilton.github.io/betaregscale/reference/logLik.brs.md)

## Examples

``` r
# Synthetic NRS-11 scores: 4 groups x 3 times. Simulated, not real data.
set.seed(2023)
nrs <- expand.grid(id = 1:80, time = c("6h", "12h", "24h"))
nrs$group <- factor(paste0("g", (nrs$id - 1) %% 4 + 1))
eta <- -1.3 + c(0, 0.75, 0.3)[nrs$time] + c(0, -0.1, 0.05, 0.1)[nrs$group]
shp <- brs_repar(mu = plogis(eta), phi = 0.3, repar = 2)
nrs$y <- round(10 * rbeta(nrow(nrs), shp$shape1, shp$shape2))

# Time only (m1) nested in time + group (m2): LR test on 3 df
m1 <- brs(y ~ time, data = nrs, ncuts = 10)
m2 <- brs(y ~ time + group, data = nrs, ncuts = 10)
anova(m1, m2)
#> Likelihood-ratio comparison of brs/brsmm models
#> 
#>          Df  logLik    AIC    BIC  Chisq Chi Df Pr(>Chisq)
#> M1 (brs)  4 -503.76 1015.5 1029.4                         
#> M2 (brs)  7 -502.64 1019.3 1043.7 2.2383      3     0.5244

# Fits under different interval directions are different response models
m2_right <- brs(y ~ time + group, data = nrs, ncuts = 10, interval = "right")
try(anova(m2, m2_right))
#> Error : All models must use the same `interval` (found: mid, right); fits with different coarsenings of the response are not comparable.
```
