# C++ log-likelihood for mixed beta interval regression (Laplace)

Computes the marginal log-likelihood for a random-intercept beta
interval model using a Laplace approximation of group-specific integrals
over Gaussian random effects.

## Usage

``` r
.brsmm_loglik_laplace_cpp(
  param,
  X,
  Z,
  y_left,
  y_right,
  yt,
  delta,
  group,
  link_mu_code,
  link_phi_code,
  repar
)
```

## Arguments

- param:

  Parameter vector: `beta (p), gamma (q), log_sigma_b (1)`.

- X:

  Mean-model design matrix.

- Z:

  Precision-model design matrix.

- y_left:

  Left interval endpoints.

- y_right:

  Right interval endpoints.

- yt:

  Midpoint responses.

- delta:

  Censoring indicator vector.

- group:

  Integer group index (1..G), one per observation.

- link_mu_code:

  Integer code for mean link.

- link_phi_code:

  Integer code for precision link.

- repar:

  Integer beta reparameterization code.

## Value

Scalar marginal log-likelihood.
