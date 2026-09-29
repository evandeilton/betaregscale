# Compute starting values for optimization

Obtains rough starting values for the beta regression parameters by
fitting a quasi-binomial GLM on the midpoint response. This provides a
reasonable initialization for the interval likelihood optimizer.

## Usage

``` r
compute_start(
  formula,
  data,
  link = NULL,
  link_phi = NULL,
  ncuts = 100L,
  lim = 0.5,
  repar = 2L,
  interval = "mid"
)
```

## Arguments

- formula:

  A [`Formula`](https://rdrr.io/pkg/Formula/man/Formula.html) object
  (possibly multi-part).

- data:

  Data frame.

- link:

  Mean link function name (`NULL`: default for `repar`).

- link_phi:

  Dispersion link function name (`NULL`: default for `repar`).

- ncuts:

  Number of scale categories.

- lim:

  Uncertainty half-width.

- repar:

  Reparameterization scheme. Under `repar = 0` the shapes are started by
  the method of moments on the midpoint response (\\f = m(1-m)/v - 1\\,
  \\p = m f\\, \\q = (1-m) f\\).

- interval:

  Interval direction used when `data` is not prepared.

## Value

Named numeric vector of starting values.
