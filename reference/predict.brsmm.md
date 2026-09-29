# Predict from a brsmm model

Predict from a brsmm model

## Usage

``` r
# S3 method for class 'brsmm'
predict(
  object,
  newdata = NULL,
  type = c("response", "link", "precision", "variance", "quantile", "score",
    "expected_score"),
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
  `"variance"`, `"quantile"`, `"score"` or `"expected_score"`. `"score"`
  is the latent score of the fit's `interval` at the conditional mean
  (support \\(0, K)\\, \\(0, K + 1)\\ or \\(-1, K)\\; about 0.5
  above/below the expected recorded score under `"right"`/`"left"`);
  `"expected_score"` is the expected recorded score \\\sum_s s\\ P(S =
  s)\\. Details:
  [`predict.brs`](https://evandeilton.github.io/betaregscale/reference/predict.brs.md).

- at:

  Numeric vector of probabilities for quantile predictions (default
  0.5).

- ...:

  Currently ignored.

## Value

Numeric vector, except when `type = "quantile"` and `at` has length
greater than 1, in which case a numeric matrix with one column per
requested quantile (named `q_<value>`, e.g. `"q_0.5"`) and one row per
observation.

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
#> [1] 0.3856376 0.3093722 0.3856376 0.3093722 0.3856376 0.5712995
head(predict(fit, type = "precision"))
#> [1] 0.3588121 0.3588121 0.3588121 0.3588121 0.3588121 0.3588121
# }
```
