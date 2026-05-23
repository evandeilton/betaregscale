# Parametric bootstrap confidence intervals for brs models

Computes bootstrap-based confidence intervals for the parameters of a
fitted `"brs"` model by repeatedly simulating data from the fitted model
and re-estimating parameters. Only `"brs"` (fixed or
variable-dispersion) objects are supported; `"brsmm"` is not supported.

## Usage

``` r
brs_bootstrap(
  object,
  R = 199L,
  level = 0.95,
  ci_type = c("percentile", "basic", "normal", "bca"),
  max_tries = NULL,
  keep_draws = FALSE
)

# S3 method for class 'brs_bootstrap'
print(x, ...)
```

## Arguments

- object:

  A fitted `"brs"` object (fixed or variable dispersion).

- R:

  Integer: number of bootstrap replicates (default 199).

- level:

  Numeric: confidence level (default 0.95).

- ci_type:

  Character: type of confidence interval. One of `"percentile"`
  (default), `"basic"`, `"normal"`, or `"bca"`.

- max_tries:

  Optional integer: maximum number of bootstrap attempts to obtain
  converged replicates. If `NULL`, uses `max(3 * R, 50)`.

- keep_draws:

  Logical: if `TRUE`, stores successful bootstrap parameter draws in
  attribute `"boot_draws"`.

- x:

  Object returned by `brs_bootstrap`.

- ...:

  Ignored.

## Value

A data frame with columns `parameter`, `estimate` (original point
estimate), `se_boot` (bootstrap standard error), `ci_lower`, `ci_upper`,
`mcse_lower`, `mcse_upper`, `wald_lower`, `wald_upper`, and `level`. The
attribute `"n_success"` gives the number of replicates that converged.
Additional attributes include `"R"`, `"n_attempted"`, `"ci_type"`, and
optionally `"boot_draws"`.

## Details

For each replicate, data are simulated via
[`brs_sim`](https://evandeilton.github.io/betaregscale/reference/brs_sim.md)
using the estimated coefficients (on the link scale) and the original
design. The model is then re-fitted with
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md).
Replicates that fail to converge are discarded; if the number of
successful replicates is too low, a warning is issued. Intervals are the
empirical quantiles of the bootstrap distribution of each parameter.

## Methods (by generic)

- `print(brs_bootstrap)`: Print method for bootstrap results

## See also

[`confint.brs`](https://evandeilton.github.io/betaregscale/reference/confint.brs.md)
for Wald intervals;
[`brs_sim`](https://evandeilton.github.io/betaregscale/reference/brs_sim.md)
for simulation;
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md) for
fitting.

## Examples

``` r
# \donttest{
dat <- data.frame(
  y = c(
    0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
    10, 40, 55, 70, 85, 25, 35, 65, 80, 15
  ),
  x1 = rep(c(1, 2), 10),
  x2 = rep(c(0, 0, 1, 1), 5)
)
prep <- brs_prep(dat, ncuts = 100)
#> brs_prep: n = 20 | exact = 0, left = 1, right = 1, interval = 18
fit <- brs(y ~ x1, data = prep)
boot <- brs_bootstrap(fit, R = 50, level = 0.95)
print(boot)
#> Bootstrap confidence intervals
#>   Level: 0.95 | CI: percentile | Successful replicates: 50 / 50 | Attempts: 50 
#> 
#>     parameter   estimate   se_boot   ci_lower   ci_upper mcse_lower mcse_upper
#> 1 (Intercept)  0.2551000 0.7772708 -0.8231021 1.93714713  0.1651500 0.15824662
#> 2          x1 -0.2202060 0.4898420 -1.3384456 0.50989159  0.1751015 0.08461005
#> 3       (phi) -0.3929144 0.3259027 -1.2359624 0.05417615  0.1536645 0.12221657
#>   wald_lower wald_upper level
#> 1 -1.4390767  1.9492767  0.95
#> 2 -1.2809286  0.8405165  0.95
#> 3 -0.9343775  0.1485488  0.95
# }
```
