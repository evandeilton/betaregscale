# Mixed Model Log-Likelihood (Eigen)

Computes marginal log-likelihood using Laplace, AGHQ, or QMC.

## Usage

``` r
brsmm_loglik_eigen(
  param,
  X,
  Z,
  Xr,
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

  \[beta, gamma, theta_re\]

- X:

  Mean design matrix

- Z:

  Precision design matrix

- Xr:

  Random-effects design matrix

- y_left, y_right, yt, delta, group:

  Data

- method:

  0=Laplace, 1=AGHQ, 2=QMC

- n_points:

  Number of quadrature/QMC points
