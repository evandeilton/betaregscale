# Log-likelihood for fixed-dispersion beta interval regression

Computes the total log-likelihood for a beta regression model with
mixed-censored responses and a single (scalar) dispersion parameter. The
heavy computation is delegated to a compiled C++ backend.

## Usage

``` r
betaregscale_loglik(
  param,
  formula,
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

The complete likelihood for observation \\i\\ with censoring indicator
\\\delta_i\\ is (Lopes, 2024, Eq. 2.24):

- \\\delta = 0\\ (uncensored):

  \\\ell_i = \log f(y_i \| a_i, b_i)\\

- \\\delta = 1\\ (left-censored):

  \\\ell_i = \log F(u_i \| a_i, b_i)\\

- \\\delta = 2\\ (right-censored):

  \\\ell_i = \log(1 - F(l_i \| a_i, b_i))\\

- \\\delta = 3\\ (interval-censored):

  \\\ell_i = \log(F(u_i \| a_i, b_i) - F(l_i \| a_i, b_i))\\

where \\a_i\\ and \\b_i\\ are the beta shape parameters derived from the
mean \\\mu_i = g^{-1}(x_i'\beta)\\ and scalar dispersion \\\phi =
h^{-1}(\gamma)\\ through the chosen reparameterization.

## Examples

``` r
set.seed(42)
n <- 100
dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
sim <- betaregscale_simulate(
  formula = ~ x1 + x2, data = dat,
  beta = c(0, 0.5, -0.2), phi = 1 / 5
)
betaregscale_loglik(
  param = c(0, 0.5, -0.2, 1 / 5),
  formula = y ~ x1 + x2, data = sim
)
#> [1] -417.4957
```
