# Predict from a fitted model

Predict from a fitted model

## Usage

``` r
# S3 method for class 'betaregscale'
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

  A fitted `"betaregscale"` object.

- newdata:

  Optional data frame for prediction.

- type:

  Prediction type: `"response"` (default), `"link"`, `"precision"`,
  `"variance"`, or `"quantile"`.

- at:

  Numeric vector of probabilities for quantile predictions (default
  0.5).

- ...:

  Ignored.

## Value

Numeric vector or matrix.
