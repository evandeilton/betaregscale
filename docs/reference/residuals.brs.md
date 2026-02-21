# Extract residuals

Extract residuals

## Usage

``` r
# S3 method for class 'brs'
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

  Currently ignored.

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

## See also

[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md),
[`fitted.brs`](https://evandeilton.github.io/betaregscale/reference/fitted.brs.md),
[`plot.brs`](https://evandeilton.github.io/betaregscale/reference/plot.brs.md)

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
head(residuals(fit))
#> [1] -0.50871087 -0.40380203 -0.30872087  0.04619797  0.24127913  0.44619797
head(residuals(fit, type = "pearson"))
#>          1          2          3          4          5          6 
#> -1.6029010 -1.2776145 -0.9727509  0.1461687  0.7602483  1.4117537 
# }
```
