# Log-likelihood for variable-dispersion beta interval regression

Computes the total log-likelihood for a beta regression model with
mixed-censored responses and observation-specific dispersion governed by
a second linear predictor. Uses the compiled C++ backend.

## Usage

``` r
betaregscale_loglik_z(
  param,
  formula = y ~ x1 + x2 | z1,
  data,
  link = "logit",
  link_phi = "logit",
  ncuts = 100L,
  type = "m",
  lim = 0.5,
  repar = 2L
)
```

## Arguments

- param:

  Numeric vector of length \\p + 1\\: the first \\p\\ elements are the
  regression coefficients \\\beta\\, and the last element is the
  (link-scale) dispersion parameter.

- formula:

  One-sided or two-sided formula for the mean model.

- data:

  Data frame containing the response and predictors.

- link:

  Character: link function for the mean (default `"logit"`).

- link_phi:

  Character: link function for the dispersion (default `"logit"`).

- ncuts:

  Integer: number of scale categories (default 100).

- type:

  **Deprecated.** Character: interval type (`"m"`, `"l"`, or `"r"`). Use
  [`bs_prepare`](https://evandeilton.github.io/betaregscale/reference/bs_prepare.md)
  to control interval geometry instead.

- lim:

  Numeric: half-width of uncertainty region (default 0.5).

- repar:

  Integer: reparameterization scheme (0, 1, or 2; default 2).

## Value

Scalar: total log-likelihood.

## Details

The formula should use the
[`Formula`](https://rdrr.io/pkg/Formula/man/Formula.html) pipe notation:
`y ~ x1 + x2 | z1 + z2`, where the left-hand side of `|` defines the
mean model and the right-hand side defines the dispersion model.

## Examples

``` r
set.seed(42)
n <- 100
dat <- data.frame(
  x1 = rnorm(n), x2 = rnorm(n),
  z1 = runif(n)
)
sim <- betaregscale_simulate_z(
  formula_x = ~ x1 + x2, formula_z = ~z1,
  data = dat,
  beta = c(0.2, -0.5, 0.3), zeta = c(0.5, -0.5)
)
betaregscale_loglik_z(
  param = c(0.2, -0.5, 0.3, 0.5, -0.5),
  formula = y ~ x1 + x2 | z1, data = sim
)
#> [1] -419.0687
```
