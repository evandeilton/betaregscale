# Simulate data from a fixed-dispersion beta interval model

Generates interval-censored observations from a beta regression model
with a single (scalar) dispersion parameter. This is useful for Monte
Carlo studies and testing.

## Usage

``` r
betaregscale_simulate(
  formula,
  data,
  beta,
  phi = 1/5,
  link = "logit",
  link_phi = "logit",
  ncuts = 100L,
  type = "m",
  lim = 0.5,
  repar = 2L
)
```

## Arguments

- formula:

  One-sided formula specifying the mean model predictors (e.g.,
  `~ x1 + x2`).

- data:

  Data frame containing the predictor variables.

- beta:

  Numeric vector of regression coefficients (length must equal the
  number of columns of the design matrix, including the intercept).

- phi:

  Scalar dispersion parameter (on the link scale).

- link:

  Mean link function (default `"logit"`).

- link_phi:

  Dispersion link function (default `"logit"`).

- ncuts:

  Number of scale categories (default 100).

- type:

  Interval type: `"m"`, `"l"`, or `"r"` (default `"m"`).

- lim:

  Half-width of the uncertainty region (default 0.5).

- repar:

  Reparameterization scheme (default 2).

## Value

A `data.frame` containing columns `left`, `right`, `yt`, `y`, and the
predictor variables.

## Examples

``` r
set.seed(42)
n <- 200
dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
sim <- betaregscale_simulate(
  formula = ~ x1 + x2, data = dat,
  beta = c(0.2, -0.5, 0.3), phi = 1 / 5,
  link = "logit", link_phi = "logit"
)
head(sim)
#>      left right      yt  y delta         x1         x2
#> 1 0.17500 0.185 0.18000 18     3  1.3709584 -2.0009292
#> 2 0.69500 0.705 0.70000 70     3 -0.5646982  0.3337772
#> 3 0.07500 0.085 0.08000  8     3  0.3631284  1.1713251
#> 4 0.61500 0.625 0.62000 62     3  0.6328626  2.0595392
#> 5 0.01500 0.025 0.02000  2     3  0.4042683 -1.3768616
#> 6 0.00001 0.005 0.00001  0     1 -0.1061245 -1.1508556
```
