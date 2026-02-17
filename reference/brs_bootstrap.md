# Parametric bootstrap confidence intervals for brs models

Computes bootstrap-based confidence intervals for the parameters of a
fitted `"brs"` model by repeatedly simulating data from the fitted model
and re-estimating parameters. Only `"brs"` (fixed or
variable-dispersion) objects are supported; `"brsmm"` is not supported.

## Usage

``` r
brs_bootstrap(object, R = 199L, level = 0.95, seed = NULL)

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

- seed:

  Optional integer: random seed for reproducibility.

- x:

  Object returned by `brs_bootstrap`.

- ...:

  Ignored.

## Value

A data frame with columns `parameter`, `estimate` (original point
estimate), `se_boot` (bootstrap standard error), `ci_lower`, `ci_upper`,
and `level`. The attribute `"n_success"` gives the number of replicates
that converged.

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

## References

Lopes, J. E. (2024). *Beta Regression for Interval-Censored
Scale-Derived Outcomes*. MSc Dissertation, PPGMNE/UFPR.

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
boot <- brs_bootstrap(fit, R = 99, level = 0.95, seed = 1)
print(boot)
#> Bootstrap confidence intervals
#>   Level: 0.95 | Successful replicates: 99 / 99 
#> 
#>     parameter   estimate   se_boot    ci_lower   ci_upper level
#> 1 (Intercept)  0.3645114 0.1614039  0.05328573  0.6213062  0.95
#> 2          x1 -0.5814425 0.1672911 -0.95176162 -0.2865863  0.95
#> 3          x2  0.4480550 0.2059549  0.08387582  0.8541168  0.95
#> 4       (phi)  0.1391969 0.1588549 -0.21002965  0.4055923  0.95
# }
```
