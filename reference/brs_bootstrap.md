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
set.seed(42)
n <- 80
dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
sim <- brs_sim(
  formula = ~ x1 + x2, data = dat,
  beta = c(0.2, -0.5, 0.3), phi = 1 / 5, ncuts = 100
)
fit <- brs(y ~ x1 + x2, data = sim)
# \donttest{
set.seed(1)  # Set seed before calling bootstrap for reproducibility
boot <- brs_bootstrap(fit, R = 99, level = 0.95)
print(boot)
#> Bootstrap confidence intervals
#>   Level: 0.95 | CI: percentile | Successful replicates: 99 / 99 | Attempts: 99 
#> 
#>     parameter   estimate   se_boot    ci_lower   ci_upper mcse_lower mcse_upper
#> 1 (Intercept)  0.3645114 0.1614039  0.05328573  0.6213062 0.04440235 0.01748219
#> 2          x1 -0.5814425 0.1672911 -0.95176162 -0.2865863 0.05662854 0.02537720
#> 3          x2  0.4480550 0.2059549  0.08387583  0.8541168 0.02833319 0.04107024
#> 4       (phi)  0.1391969 0.1588549 -0.21002965  0.4055923 0.08952651 0.04770865
#>    wald_lower wald_upper level
#> 1  0.05918412  0.6698387  0.95
#> 2 -0.88063079 -0.2822543  0.95
#> 3  0.09586929  0.8002408  0.95
#> 4 -0.13974933  0.4181431  0.95
# }
```
