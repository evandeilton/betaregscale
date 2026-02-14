# Fit a variable-dispersion beta interval regression model

Estimates the parameters of a beta regression model with
observation-specific dispersion governed by a second linear predictor.
Both submodels are estimated jointly via maximum likelihood, using the
complete likelihood with mixed censoring (Lopes, 2024, Eq. 2.24).

## Usage

``` r
betaregscale_fit_z(
  formula,
  data,
  link = "logit",
  link_phi = "logit",
  hessian_method = c("numDeriv", "optim"),
  ncuts = 100L,
  type = "m",
  lim = 0.5,
  repar = 2L,
  method = c("BFGS", "L-BFGS-B")
)
```

## Arguments

- formula:

  A [`Formula`](https://rdrr.io/pkg/Formula/man/Formula.html)-style
  formula with two parts: `y ~ x1 + x2 | z1 + z2`.

- data:

  Data frame.

- link:

  Mean link function (default `"logit"`).

- link_phi:

  Dispersion link function (default `"logit"`).

- hessian_method:

  Character: `"numDeriv"` or `"optim"`.

- ncuts:

  Number of scale categories (default 100).

- type:

  **Deprecated.** Interval type (default `"m"`). Use
  [`bs_prepare`](https://evandeilton.github.io/betaregscale/reference/bs_prepare.md)
  to control interval geometry instead.

- lim:

  Uncertainty half-width (default 0.5).

- repar:

  Reparameterization scheme (default 2).

- method:

  Optimization method (default `"BFGS"`).

## Value

An object of class `"betaregscale"`.

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
  beta = c(0.2, -0.5, 0.3),
  zeta = c(1, 1.2)
)
fit <- betaregscale_fit_z(
  formula = y ~ x1 + x2 | z1, data = sim,
  link = "logit", link_phi = "logit"
)
print(fit)
#> 
#> Call:
#> betaregscale_fit_z(formula = y ~ x1 + x2 | z1, data = sim, link = "logit", 
#>     link_phi = "logit")
#> 
#> Coefficients (mean model with logit link):
#> (Intercept)          x1          x2 
#>      0.0959     -0.3912      0.2965 
#> 
#> Phi coefficients (precision model with logit link):
#> (phi)_(Intercept)          (phi)_z1 
#>            1.4876            0.0040 
#> 
```
