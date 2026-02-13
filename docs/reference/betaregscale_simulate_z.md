# Simulate data from a variable-dispersion beta interval model

Generates interval-censored observations from a beta regression model
with observation-specific dispersion governed by a second linear
predictor.

## Usage

``` r
betaregscale_simulate_z(
  formula_x = ~x1 + x2,
  formula_z = ~z1 + z2,
  data,
  beta = c(0, 0.5, -0.2),
  zeta = c(1, 0.5, 0.2),
  link = "logit",
  link_phi = "logit",
  ncuts = 100L,
  type = "m",
  lim = 0.5,
  repar = 2L
)
```

## Arguments

- formula_x:

  One-sided formula for the mean model predictors.

- formula_z:

  One-sided formula for the dispersion model predictors.

- data:

  Data frame containing the predictor variables.

- beta:

  Numeric vector of mean-model coefficients.

- zeta:

  Numeric vector of dispersion-model coefficients.

- link:

  Mean link function (default `"logit"`).

- link_phi:

  Dispersion link function (default `"logit"`).

- ncuts:

  Number of scale categories (default 100).

- type:

  Interval type (default `"m"`).

- lim:

  Uncertainty half-width (default 0.5).

- repar:

  Reparameterization scheme (default 2).

## Value

A `data.frame` with interval endpoints and predictors.

## Examples

``` r
set.seed(42)
n <- 200
dat <- data.frame(
  x1 = rnorm(n), x2 = rnorm(n),
  z1 = runif(n), z2 = runif(n)
)
sim <- betaregscale_simulate_z(
  formula_x = ~ x1 + x2, formula_z = ~ z1 + z2,
  data = dat,
  beta = c(0.2, -0.5, 0.3),
  zeta = c(0.2, -0.4, 0.2)
)
head(sim)
#>    left right   yt  y delta         x1         x2        z1         z2
#> 1 0.615 0.625 0.62 62     3  1.3709584 -2.0009292 0.9090475 0.84829322
#> 2 0.895 0.905 0.90 90     3 -0.5646982  0.3337772 0.8999248 0.06274633
#> 3 0.985 0.995 0.99 99     3  0.3631284  1.1713251 0.1923493 0.81984509
#> 4 0.985 0.995 0.99 99     3  0.6328626  2.0595392 0.5322903 0.53936029
#> 5 0.865 0.875 0.87 87     3  0.4042683 -1.3768616 0.5221247 0.49902010
#> 6 0.255 0.265 0.26 26     3 -0.1061245 -1.1508556 0.1603357 0.02222732
```
