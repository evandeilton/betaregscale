# Mixed Model Log-Likelihood (Eigen)

Computes marginal log-likelihood using Laplace, AGHQ, or QMC.

## Usage

``` r
brsmm_loglik_eigen(
  param,
  X,
  Z,
  y_left,
  y_right,
  yt,
  delta,
  group,
  link_mu,
  link_phi,
  repar,
  method,
  n_points
)
```

## Arguments

- param:

  \[beta, gamma, log_sigma\]

- X, :

  Z Design matrices

- y_left, :

  y_right, yt, delta Data

- group:

  Group indices

- method:

  0=Laplace, 1=AGHQ, 2=QMC

- n_points:

  Number of quadrature/QMC points
