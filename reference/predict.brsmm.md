# Predict from a brsmm model

Predict from a brsmm model

## Usage

``` r
# S3 method for class 'brsmm'
predict(
  object,
  newdata = NULL,
  type = c("response", "link", "precision", "variance"),
  ...
)
```

## Arguments

- object:

  A fitted `"brsmm"` object.

- newdata:

  Optional data frame.

- type:

  Character: `"response"` (default), `"link"`, `"precision"`, or
  `"variance"`.

- ...:

  Currently ignored.

## Value

Numeric vector.
