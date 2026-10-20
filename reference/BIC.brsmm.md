# BIC for brsmm models

\\-2\ell + \log(n)\\k\\, with \\k\\ the number of parameters (fixed
effects, precision, packed random-effect and cutpoint parameters) and
\\n\\ the number of observations (`n = "obs"`, the default, as `lme4`
does) or of groups (`n = "groups"`). There is no single sample size for
a mixed model; counting groups penalizes less when groups have several
observations each. Compare BIC values only between fits with the same
convention.

## Usage

``` r
# S3 method for class 'brsmm'
BIC(object, ..., n = c("obs", "groups"))
```

## Arguments

- object:

  A fitted `"brsmm"` object.

- ...:

  Further fitted models; with `n = "groups"` they must be `"brsmm"`
  fits, otherwise their own [`BIC()`](https://rdrr.io/r/stats/AIC.html)
  is used.

- n:

  `"obs"` (default; `nobs(object)`) or `"groups"` (`object$ngroups`).

## Value

Numeric scalar; with several models, a data frame with columns `df` and
`BIC`, one row per model, as
[`stats::BIC()`](https://rdrr.io/r/stats/AIC.html).

## See also

[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md),
[`logLik.brsmm`](https://evandeilton.github.io/betaregscale/reference/logLik.brsmm.md),
[`AIC.brsmm`](https://evandeilton.github.io/betaregscale/reference/AIC.brsmm.md),
[`brs_gof`](https://evandeilton.github.io/betaregscale/reference/brs_gof.md)

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
BIC(fit)
#> [1] 196.3492
BIC(fit, n = "groups")
#> [1] 189.9115
# }
```
