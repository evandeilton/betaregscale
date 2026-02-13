# Extract model coefficients

Extract model coefficients

## Usage

``` r
# S3 method for class 'betaregscale'
coef(object, model = c("full", "mean", "precision"), ...)
```

## Arguments

- object:

  A fitted `"betaregscale"` object.

- model:

  Character: which component to return. `"full"` (default) returns all
  parameters, `"mean"` returns only the mean-model coefficients,
  `"precision"` returns only the precision coefficients.

- ...:

  Ignored.

## Value

Named numeric vector of estimated parameters.
