# C++ gradient for fixed-dispersion log-likelihood

Returns the gradient vector of the log-likelihood with respect to all
parameters (beta coefficients + scalar phi), using a central-difference
numerical approximation (step = 1e-6).

## Usage

``` r
.betaregscale_grad_fixed_cpp(
  param,
  X,
  y_left,
  y_right,
  yt,
  delta,
  link_mu_code,
  link_phi_code,
  repar
)
```

## Arguments

- param:

  Parameter vector (same layout as loglik function).

- X:

  Design matrix (n x p).

- y_left:

  Left endpoints.

- y_right:

  Right endpoints.

- yt:

  Midpoint responses.

- delta:

  Integer censoring indicators.

- link_mu_code:

  Integer mean link code.

- link_phi_code:

  Integer dispersion link code.

- repar:

  Integer reparameterization type.

## Value

Numeric gradient vector of length `ncol(X) + 1`.
