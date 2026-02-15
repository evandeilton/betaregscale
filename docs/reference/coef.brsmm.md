# Extract coefficients from a brsmm fit

Extract coefficients from a brsmm fit

## Usage

``` r
# S3 method for class 'brsmm'
coef(object, model = c("full", "mean", "precision", "random"), ...)
```

## Arguments

- object:

  A fitted `"brsmm"` object.

- model:

  Character: `"full"` (default), `"mean"`, `"precision"`, or `"random"`.

- ...:

  Ignored.

## Value

Named numeric vector.
