# Summarize a fitted brsmm model

Wald tests for the fixed effects. The random effects are reported as
standard deviations and correlations (`varcorr`) with Wald intervals
built on a transformed scale and mapped back: \\\exp\\ of the interval
for \\\log SD\\, \\\tanh\\ of the interval for \\\mathrm{atanh}(\rho)\\
(delta method from the packed Cholesky parameters). No test or p-value
is given for them: a z-test of \\\log SD\\ tests \\SD = 1\\, and \\SD =
0\\ lies on the boundary; use
[`anova.brsmm`](https://evandeilton.github.io/betaregscale/reference/anova.brsmm.md)
(chi-bar-square mixture) against the model without the term. The
randomized quantile residuals are drawn without changing the caller's
RNG state.

## Usage

``` r
# S3 method for class 'brsmm'
summary(object, level = 0.95, ...)
```

## Arguments

- object:

  A fitted `"brsmm"` object.

- level:

  Confidence level of the `varcorr` intervals.

- ...:

  Currently ignored.

## Value

Object of class `"summary.brsmm"`; `coefficients$random` holds the
packed Cholesky parameters (estimate and standard error only) and
`varcorr` the SD/correlation table.

## See also

[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md),
[`print.summary.brsmm`](https://evandeilton.github.io/betaregscale/reference/print.summary.brsmm.md),
[`brs_gof`](https://evandeilton.github.io/betaregscale/reference/brs_gof.md),
[`brsmm_re_study`](https://evandeilton.github.io/betaregscale/reference/brsmm_re_study.md)

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
s <- summary(fit)
s$coefficients$mean
#>               Estimate Std. Error    z value  Pr(>|z|)
#> (Intercept)  0.4211917  0.7770316  0.5420522 0.5877826
#> x1          -0.3373687  0.4752862 -0.7098222 0.4778144
# }
```
