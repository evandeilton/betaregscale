# Predict from a brsmm model

Predict from a brsmm model

## Usage

``` r
# S3 method for class 'brsmm'
predict(
  object,
  newdata = NULL,
  type = c("response", "link", "precision", "variance", "quantile"),
  at = 0.5,
  ...
)
```

## Arguments

- object:

  A fitted `"brsmm"` object.

- newdata:

  Optional data frame.

- type:

  Character: `"response"` (default), `"link"`, `"precision"`,
  `"variance"`, or `"quantile"`.

- at:

  Numeric vector of probabilities for quantile predictions (default
  0.5).

- ...:

  Currently ignored.

## Value

Numeric vector.

## See also

[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md),
[`fitted.brsmm`](https://evandeilton.github.io/betaregscale/reference/fitted.brsmm.md),
[`brs_predict_scoreprob`](https://evandeilton.github.io/betaregscale/reference/brs_predict_scoreprob.md)

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
head(predict(fit))
#>         1         2         3         4         5         6 
#> 0.3857010 0.3093579 0.3857010 0.3093579 0.3857010 0.5712864 
head(predict(fit, type = "precision"))
#> [1] 0.3588217 0.3588217 0.3588217 0.3588217 0.3588217 0.3588217
# }
```
