# Extract random effects from a brsmm model

Extract random effects from a brsmm model

## Usage

``` r
# S3 method for class 'brsmm'
ranef(object, ...)
```

## Arguments

- object:

  A fitted `"brsmm"` object.

- ...:

  Currently ignored.

## Value

A matrix or named numeric vector of group-specific random-effect
posterior modes.

## See also

[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md),
[`brsmm_re_study`](https://evandeilton.github.io/betaregscale/reference/brsmm_re_study.md),
[`ranef`](https://evandeilton.github.io/betaregscale/reference/ranef.md)

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
ranef(fit)
#>           1           2           3           4 
#> -0.54955383  0.54074901  0.03720833 -0.02380723 
# }
```
