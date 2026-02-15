# Simulate data from beta interval models

Simulates interval-censored responses from fixed- or variable-dispersion
beta regression models.

## Usage

``` r
brs_sim(
  formula,
  data,
  beta,
  phi = 1/5,
  zeta = NULL,
  link = "logit",
  link_phi = "logit",
  ncuts = 100L,
  lim = 0.5,
  repar = 2L,
  delta = NULL
)
```

## Arguments

- formula:

  Model formula with one (mean) or two parts (mean `|` precision). A
  left-hand-side response is allowed but ignored.

- data:

  Data frame with predictor variables.

- beta:

  Numeric vector of mean-model coefficients.

- phi:

  Scalar dispersion parameter (link scale), used only for one-part
  formulas.

- zeta:

  Numeric vector of precision-model coefficients (link scale), required
  for two-part formulas.

- link:

  Mean link function.

- link_phi:

  Precision link function.

- ncuts:

  Number of scale categories.

- lim:

  Half-width used in interval construction.

- repar:

  Reparameterization scheme.

- delta:

  Forced censoring type (`0,1,2,3`) or `NULL`.

## Value

A data frame with columns `left`, `right`, `yt`, `y`, `delta`, plus
simulated predictor columns from the model matrices. When
`delta != NULL`, the output carries `attr(, "is_prepared") = TRUE`.

## Details

The model structure is controlled by `formula` in the same style as
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md):

- one-part formula (`~ x1 + x2` or `y ~ x1 + x2`): fixed dispersion
  using scalar `phi`.

- two-part formula (`~ x1 + x2 | z1` or `y ~ x1 + x2 | z1`): variable
  dispersion using coefficient vector `zeta`.

The `delta` argument can force a single censoring type (`0,1,2,3`) for
all observations; otherwise, censoring is classified automatically from
simulated scale values via
[`brs_check`](https://evandeilton.github.io/betaregscale/reference/brs_check.md).

## Examples

``` r
set.seed(42)
n <- 200
dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n), z1 = runif(n))

# Fixed dispersion
sim_fixed <- brs_sim(
  formula = ~ x1 + x2, data = dat,
  beta = c(0.2, -0.5, 0.3), phi = 1 / 5
)

# Variable dispersion
sim_var <- brs_sim(
  formula = ~ x1 + x2 | z1, data = dat,
  beta = c(0.2, -0.5, 0.3), zeta = c(0.5, -0.5)
)
```
