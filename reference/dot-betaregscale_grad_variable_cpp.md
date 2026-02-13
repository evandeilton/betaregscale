# C++ gradient for variable-dispersion log-likelihood

Central-difference gradient for the variable-dispersion model.

## Usage

``` r
.betaregscale_grad_variable_cpp(
  param,
  X,
  Z,
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

  Parameter vector.

- X:

  Mean design matrix.

- Z:

  Dispersion design matrix.

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

Numeric gradient vector.
