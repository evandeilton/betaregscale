# Fit a fixed-dispersion beta interval regression model

Estimates the parameters of a beta regression model with a single
(scalar) dispersion parameter using maximum likelihood. The
log-likelihood and its gradient are evaluated by the compiled C++
backend supporting the complete likelihood with mixed censoring types
(Lopes, 2024, Eq. 2.24).

## Usage

``` r
brs_fit_fixed(
  formula,
  data,
  link = "logit",
  link_phi = "logit",
  ncuts = 100L,
  lim = 0.5,
  hessian_method = c("numDeriv", "optim"),
  repar = 2L,
  method = c("BFGS", "L-BFGS-B")
)
```

## Arguments

- formula:

  Two-sided formula `y ~ x1 + x2 + ...`.

- data:

  Data frame.

- link:

  Mean link function (default `"logit"`).

- link_phi:

  Dispersion link function (default `"logit"`).

- ncuts:

  Number of scale categories (default 100).

- lim:

  Uncertainty half-width (default 0.5).

- hessian_method:

  Character: `"numDeriv"` (default) or `"optim"`. With `"numDeriv"` the
  Hessian is computed after convergence using
  [`hessian`](https://rdrr.io/pkg/numDeriv/man/hessian.html), which is
  typically more accurate than the built-in optim Hessian.

- repar:

  Reparameterization scheme (default 2).

- method:

  Optimization method: `"BFGS"` (default) or `"L-BFGS-B"`.

## Value

An object of class `"brs"`.

## Examples

``` r
set.seed(42)
n <- 100
dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
sim <- brs_sim(
  formula = ~ x1 + x2, data = dat,
  beta = c(0.2, -0.5, 0.3), phi = 1 / 5
)
fit <- brs_fit_fixed(
  formula = y ~ x1 + x2, data = sim,
  link = "logit", link_phi = "logit"
)
print(fit)
#> 
#> Call:
#> brs_fit_fixed(formula = y ~ x1 + x2, data = sim, link = "logit", 
#>     link_phi = "logit")
#> 
#> Coefficients (mean model with logit link):
#> (Intercept)          x1          x2 
#>      0.0969     -0.5117      0.1147 
#> 
#> Phi coefficients (precision model with logit link):
#>  (phi) 
#> 0.1612 
#> 
```
