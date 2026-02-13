# Wald confidence intervals

Computes Wald confidence intervals for model parameters using the normal
approximation (Lopes, 2024, Eq. 2.30–2.31).

## Usage

``` r
# S3 method for class 'betaregscale'
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

  Ignored.

## Value

Matrix with columns for lower and upper confidence bounds.
