# Simulate data from a fixed-dispersion beta interval model

Generates observations from a beta regression model with a single
(scalar) dispersion parameter. This is useful for Monte Carlo studies,
power analysis, and testing. The `delta` argument controls the censoring
type of the simulated data.

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
  repar = 2L,
  delta = NULL
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

  Mean link function (default `"logit"`). Supported: `"logit"`,
  `"probit"`, `"cloglog"`, `"cauchit"`, `"log"`.

- link_phi:

  Dispersion link function (default `"logit"`).

- ncuts:

  Integer: number of scale categories \\K\\ (default 100).

- type:

  **Deprecated.** Interval type: `"m"`, `"l"`, or `"r"`. Use
  [`bs_prepare`](https://evandeilton.github.io/betaregscale/reference/bs_prepare.md)
  to control interval geometry instead.

- lim:

  Numeric: half-width \\h\\ of the uncertainty region (default 0.5).

- repar:

  Integer: reparameterization scheme (default 2). 0 = direct \\(a, b)\\,
  1 = precision/Ferrari \\(\mu, \phi = a + b)\\, 2 = mean-variance/Bayer
  \\(\mu, \phi = 1/(a + b + 1))\\.

- delta:

  Integer or `NULL`. If `NULL` (default), censoring is determined
  automatically by
  [`check_response`](https://evandeilton.github.io/betaregscale/reference/check_response.md).

  If an integer in `{0, 1, 2, 3}`, **all** simulated observations are
  forced to that censoring type. The actual simulated values are
  preserved so that observation-specific endpoints reflect the
  underlying covariate-driven variation. See Details.

## Value

A `data.frame` with \\n\\ rows and columns: `left`, `right`, `yt`, `y`,
`delta`, plus the predictor columns from `data`. When `delta != NULL`,
the data frame carries the attribute `"bs_prepared" = TRUE`.

## Details

**Data generation process**:

1.  The design matrix \\X\\ is built from `formula` and `data`.

2.  The linear predictor is \\\eta = X \beta\\, and the mean is \\\mu =
    g^{-1}(\eta)\\ where \\g\\ is the link function.

3.  The scalar dispersion is \\\phi = h^{-1}(\texttt{phi})\\ where \\h\\
    is the dispersion link.

4.  Beta shape parameters \\(a, b)\\ are derived from \\(\mu, \phi)\\
    via the chosen reparameterization scheme.

5.  Raw values \\y^\*\_i \sim \text{Beta}(a_i, b_i)\\ are drawn on \\(0,
    1)\\.

6.  The raw values are transformed into the response matrix by
    `.build_simulated_response()` (see below).

**Role of the `delta` argument**:

When `delta = NULL` (default), the raw values are rounded to the scale
grid (\\y\_{\text{grid}} = \text{round}(y^\* \times K)\\) and passed to
[`check_response`](https://evandeilton.github.io/betaregscale/reference/check_response.md)
for automatic classification: \\y = 0 \to \delta = 1\\, \\y = K \to
\delta = 2\\, otherwise \\\delta = 3\\. The resulting dataset has a
natural mix of censoring types driven by the simulated values.

When `delta` is an integer in \\\\0, 1, 2, 3\\\\, **all** observations
are forced to that censoring type, but the actual simulated \\y^\*\\
values are **preserved** on the grid so that each observation retains
its covariate-driven variation. Specifically:

- `delta = 0` (exact):

  The continuous \\y^\*\\ values are used directly on \\(0, 1)\\; \\l_i
  = u_i = y_t = y^\*\_i\\.

- `delta = 1` (left-censored):

  The grid values \\y\_{\text{grid}} = \text{round}(y^\* K)\\ are kept.
  [`check_response()`](https://evandeilton.github.io/betaregscale/reference/check_response.md)
  is called with forced `delta = rep(1, n)`, producing: \\l_i =
  \epsilon\\, \\u_i = (y\_{\text{grid}} + h) / K\\ for non-boundary, or
  \\u_i = h / K\\ when \\y\_{\text{grid}} = 0\\.

- `delta = 2` (right-censored):

  Same logic: \\u_i = 1 - \epsilon\\, \\l_i = (y\_{\text{grid}} - h) /
  K\\ for non-boundary, or \\l_i = (K - h) / K\\ when \\y\_{\text{grid}}
  = K\\.

- `delta = 3` (interval-censored):

  Grid values are clamped to \\\[1, K-1\]\\ (avoiding boundaries) and
  [`check_response()`](https://evandeilton.github.io/betaregscale/reference/check_response.md)
  is called with forced `delta = rep(3, n)`.

**Attribute `"bs_prepared"`**:

When `delta != NULL`, the returned data frame carries the attribute
`"bs_prepared" = TRUE`. This signals to `.extract_response()` (and thus
to
[`betaregscale`](https://evandeilton.github.io/betaregscale/reference/betaregscale.md),
[`betaregscale_loglik`](https://evandeilton.github.io/betaregscale/reference/betaregscale_loglik.md),
etc.) that the pre-computed columns `left`, `right`, `yt`, and `delta`
should be used directly, bypassing the automatic classification of
[`check_response`](https://evandeilton.github.io/betaregscale/reference/check_response.md).
Without this attribute, the fitting functions would re-classify the
response from the `y` column alone, which would ignore the forced delta
and produce incorrect censoring indicators (e.g., an observation with
\\y = 50\\ and forced \\\delta = 2\\ would be reclassified as \\\delta =
3\\ by the boundary rules).

When `delta = NULL`, the attribute is **not** set, so the default
pipeline applies.

## See also

[`betaregscale_simulate_z`](https://evandeilton.github.io/betaregscale/reference/betaregscale_simulate_z.md)
for variable- dispersion simulation;
[`check_response`](https://evandeilton.github.io/betaregscale/reference/check_response.md)
for the endpoint computation rules;
[`bs_prepare`](https://evandeilton.github.io/betaregscale/reference/bs_prepare.md)
for analyst-facing data pre-processing.

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

# Force all observations to be interval-censored
sim3 <- betaregscale_simulate(
  formula = ~ x1 + x2, data = dat,
  beta = c(0.2, -0.5, 0.3), phi = 1 / 5,
  delta = 3
)
table(sim3$delta)
#> 
#>   3 
#> 200 

# Force right-censored: y values vary, all delta = 2
sim2 <- betaregscale_simulate(
  formula = ~ x1 + x2, data = dat,
  beta = c(0.2, -0.5, 0.3), phi = 1 / 5,
  delta = 2
)
head(sim2[, c("left", "right", "y", "delta")])
#>    left   right   y delta
#> 1 0.005 0.99999   1     2
#> 2 0.495 0.99999  50     2
#> 3 0.015 0.99999   2     2
#> 4 0.995 0.99999 100     2
#> 5 0.015 0.99999   2     2
#> 6 0.115 0.99999  12     2
# Note: left varies per observation, right = 1 - eps
```
