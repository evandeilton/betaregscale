# Wald confidence intervals

Computes Wald confidence intervals for model parameters using the normal
approximation.

## Usage

``` r
# S3 method for class 'brs'
confint(
  object,
  parm,
  level = 0.95,
  model = c("full", "mean", "precision"),
  ...
)
```

## Arguments

- object:

  A fitted `"betaregscale"` object.

- parm:

  Character or integer: which parameters. If missing, all parameters are
  returned.

- level:

  Confidence level (default 0.95).

- model:

  Character: `"full"`, `"mean"`, or `"precision"`.

- ...:

  Currently ignored.

## Value

Matrix with columns for lower and upper confidence bounds.
