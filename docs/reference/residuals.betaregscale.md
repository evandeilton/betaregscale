# Extract residuals

Extract residuals

## Usage

``` r
# S3 method for class 'betaregscale'
residuals(
  object,
  type = c("response", "pearson", "deviance", "rqr", "weighted", "sweighted"),
  ...
)
```

## Arguments

- object:

  A fitted `"betaregscale"` object.

- type:

  Residual type. One of `"response"` (default), `"pearson"`,
  `"deviance"`, `"rqr"` (randomized quantile), `"weighted"`, or
  `"sweighted"`.

- ...:

  Ignored.

## Value

Numeric vector of residuals.

## Details

For Pearson residuals the variance formula depends on the
reparameterization stored in `object$repar`:

- repar = 1 (precision):

  V = mu(1 - mu) / (1 + phi)

- repar = 2 (mean-variance):

  V = mu(1 - mu) \* phi

The weighted and sweighted residuals use the digamma/trigamma
formulation from the precision parameterization (repar = 1), so internal
conversion is applied when `repar != 1`.
