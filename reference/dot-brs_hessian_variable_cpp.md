# C++ Hessian for the variable-dispersion log-likelihood

Chain rule on the linear predictors: blocks `crossprod(X, w_mm * X)`,
`crossprod(X, w_mp * Z)` and `crossprod(Z, w_pp * Z)`.

## Usage

``` r
.brs_hessian_variable_cpp(
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

  Numeric vector: `ncol(X)` beta then `ncol(Z)` gamma.

- X, Z:

  Design matrices of the mean (n x p) and dispersion (n x q).

- y_left, y_right:

  Interval endpoints on (0, 1).

- yt:

  Exact response on (0, 1) (used when `delta = 0`).

- delta:

  Integer censoring indicators (0,1,2,3).

- link_mu_code, link_phi_code:

  Integer link codes.

- repar:

  Integer reparameterization type (0, 1, or 2).

## Value

Symmetric matrix of order `ncol(X) + ncol(Z)` (log-likelihood scale).
