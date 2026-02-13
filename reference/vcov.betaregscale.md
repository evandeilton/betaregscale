# Variance-covariance matrix of estimated coefficients

Variance-covariance matrix of estimated coefficients

## Usage

``` r
# S3 method for class 'betaregscale'
vcov(object, model = c("full", "mean", "precision"), ...)
```

## Arguments

- object:

  A fitted `"betaregscale"` object.

- model:

  Character: which component (`"full"`, `"mean"`, or `"precision"`).

- ...:

  Ignored.

## Value

A square numeric matrix.
