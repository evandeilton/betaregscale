# Apply the inverse-link function to a linear predictor

Evaluates the inverse of a standard link function for a given
linear-predictor vector or scalar. This is a convenience wrapper around
[`make.link`](https://rdrr.io/r/stats/make.link.html).

## Usage

``` r
apply_inv_link(eta, link)
```

## Arguments

- eta:

  Numeric vector or scalar — the linear predictor \\\eta = X \beta\\.

- link:

  Character string naming the link function. Supported values:
  `"logit"`, `"probit"`, `"cauchit"`, `"cloglog"`, `"log"`, `"sqrt"`,
  `"1/mu^2"`, `"inverse"`, `"identity"`.

## Value

Numeric vector (or scalar) of the same length as `eta`, containing
\\g^{-1}(\eta)\\.
