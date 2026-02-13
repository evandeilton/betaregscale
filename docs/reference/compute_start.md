# Compute starting values for optimization

Obtains rough starting values for the beta regression parameters by
fitting a quasi-binomial GLM on the midpoint response. This provides a
reasonable initialization for the interval likelihood optimizer.

## Usage

``` r
compute_start(
  formula,
  data,
  link = "logit",
  link_phi = "logit",
  ncuts = 100L,
  type = "m",
  lim = 0.5
)
```

## Arguments

- formula:

  A [`Formula`](https://rdrr.io/pkg/Formula/man/Formula.html) object
  (possibly multi-part).

- data:

  Data frame.

- link:

  Mean link function name.

- link_phi:

  Dispersion link function name.

- ncuts:

  Number of scale categories.

- type:

  Interval type.

- lim:

  Uncertainty half-width.

## Value

Named numeric vector of starting values.
