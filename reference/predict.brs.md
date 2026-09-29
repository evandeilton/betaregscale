# Predict from a fitted model

Predict from a fitted model

## Usage

``` r
# S3 method for class 'brs'
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

  A fitted `"brs"` object.

- newdata:

  Optional data frame for prediction.

- type:

  Prediction type: `"response"` (default; the mean \\E\[Y\] = a / (a +
  b)\\), `"link"` (linear predictor of the first parameter),
  `"precision"` (second parameter on its own scale), `"variance"`, or
  `"quantile"`.

- at:

  Numeric vector of probabilities for quantile predictions (default
  0.5).

- ...:

  Currently ignored.

## Value

Numeric vector or matrix.

## See also

[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md),
[`fitted.brs`](https://evandeilton.github.io/betaregscale/reference/fitted.brs.md),
[`brs_predict_scoreprob`](https://evandeilton.github.io/betaregscale/reference/brs_predict_scoreprob.md)

## Examples

``` r
# \donttest{
dat <- data.frame(
  y = c(
    0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
    10, 40, 55, 70, 85, 25, 35, 65, 80, 15
  ),
  x1 = rep(c(1, 2), 10)
)
prep <- brs_prep(dat, ncuts = 100)
#> brs_prep: n = 20 | exact = 0, left = 1, right = 1, interval = 18
fit <- brs(y ~ x1, data = prep)
head(predict(fit))
#> [1] 0.5087226 0.4538041 0.5087226 0.4538041 0.5087226 0.4538041
head(predict(fit, type = "precision"))
#> [1] 0.4030159 0.4030159 0.4030159 0.4030159 0.4030159 0.4030159
newdat <- data.frame(x1 = c(1, 2))
predict(fit, newdata = newdat)
#> [1] 0.5087226 0.4538041
# }
```
