# C++ gradient for the fixed-dispersion log-likelihood

Chain rule on the linear predictors: `crossprod(X, d_mu)` and
`sum(d_phi)`, with per-observation derivatives by Richardson central
differences (brs_deriv.h).

## Usage

``` r
.brs_grad_fixed_cpp(
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

  Numeric vector: `ncol(X)` coefficients, then phi (link scale).

- X:

  Design matrix (n x p).

- y_left, y_right:

  Interval endpoints on (0, 1).

- yt:

  Exact response on (0, 1) (used when `delta = 0`).

- delta:

  Integer censoring indicators (0,1,2,3).

- link_mu_code, link_phi_code:

  Integer link codes (see `link_to_code`).

- repar:

  Integer reparameterization type (0, 1, or 2).

## Value

Numeric gradient vector of length `ncol(X) + 1`.
