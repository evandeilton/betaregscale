# Fit a beta interval regression model

Unified interface that dispatches to
[`betaregscale_fit`](https://evandeilton.github.io/betaregscale/reference/betaregscale_fit.md)
(fixed dispersion) or
[`betaregscale_fit_z`](https://evandeilton.github.io/betaregscale/reference/betaregscale_fit_z.md)
(variable dispersion) based on the formula structure.

## Usage

``` r
betaregscale(
  formula,
  data,
  link = "logit",
  link_phi = "logit",
  ncuts = 100L,
  type = "m",
  lim = 0.5,
  repar = 2L,
  method = c("BFGS", "L-BFGS-B"),
  hessian_method = c("numDeriv", "optim")
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

- hessian_method:

  Character: `"numDeriv"` or `"optim"`.

## Value

An object of class `"betaregscale"`.

## Details

If the formula contains a `|` separator (e.g., `y ~ x1 + x2 | z1`), the
variable-dispersion model is fitted; otherwise, a fixed-dispersion model
is used.

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

# Fixed dispersion
fit1 <- betaregscale(y ~ x1 + x2, data = sim)
#> Warning: The 'type' argument of betaregscale_fit() is deprecated and will be removed in a future version. Use bs_prepare() to control interval geometry.
print(fit1)
#> 
#> Call:
#> betaregscale(formula = y ~ x1 + x2, data = sim)
#> 
#> Coefficients (mean model with logit link):
#> (Intercept)          x1          x2 
#>      0.0959     -0.3912      0.2964 
#> 
#> Phi coefficients (precision model with logit link):
#>  (phi) 
#> 1.4895 
#> 

# Variable dispersion
fit2 <- betaregscale(y ~ x1 + x2 | z1, data = sim)
#> Warning: The 'type' argument of betaregscale_fit_z() is deprecated and will be removed in a future version. Use bs_prepare() to control interval geometry.
print(fit2)
#> 
#> Call:
#> betaregscale(formula = y ~ x1 + x2 | z1, data = sim)
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
