# Apply the forward link function to a response value

Evaluates the link function \\g(\mu)\\ for a given response vector. Used
internally for starting-value computation.

## Usage

``` r
apply_link(mu, link)
```

## Arguments

- mu:

  Numeric vector of response values.

- link:

  Character string naming the link function (same set as
  [`apply_inv_link`](https://evandeilton.github.io/betaregscale/reference/apply_inv_link.md)).

## Value

Numeric vector containing \\g(\mu)\\.
