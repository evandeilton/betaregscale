# Group-level modes for random intercept (Laplace)

Computes the posterior mode and local SD for each group random intercept
under the current parameter vector.

## Usage

``` r
.brsmm_group_modes_cpp(
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

Numeric matrix with columns `mode_b` and `sd_b`.
