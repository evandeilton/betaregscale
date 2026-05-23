# Summarize a fitted brsmm model

Summarize a fitted brsmm model

## Usage

``` r
# S3 method for class 'brsmm'
summary(object, ...)
```

## Arguments

- object:

  A fitted `"brsmm"` object.

- ...:

  Currently ignored.

## Value

Object of class `"summary.brsmm"`.

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
#> (Intercept)  0.4213505  0.9436722  0.4465009 0.6552355
#> x1          -0.3374768  0.5737245 -0.5882210 0.5563840
# }
```
