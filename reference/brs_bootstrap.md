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
  (default), `"basic"`, `"normal"`, or `"bca"`. See the section on the
  cost of `"bca"` below.

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
Additional attributes include `"R"`, `"n_attempted"`, `"n_failed"`,
`"fail_rate"`, `"n_jack_failed"` (BCa only), `"ci_type"`, and optionally
`"boot_draws"`.

## Details

Each replicate draws a new response \\y^\*\_i \sim \mathrm{Beta}(a_i,
b_i)\\ at the fitted shapes of the rows used by the fit and substitutes
it into a copy of `object$data`; covariates, factor levels,
transformations and the formula stay those of the original fit, which is
then re-fitted with
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md)
under the same `ncuts`, `lim`, `interval`, links, `repar` and optimizer.
The observation mechanism of each row is reproduced: exact observations
(\\\delta = 0\\) stay continuous; scores are re-coarsened on the fit's
grid (the mapping of
[`brs_sim`](https://evandeilton.github.io/betaregscale/reference/brs_sim.md)).
Rows whose bounds came from the analyst
([`brs_prep`](https://evandeilton.github.io/betaregscale/reference/brs_prep.md)
Modes 2–4, bounds that are not the cell of the row's score) keep their
thresholds as fixed and non-informative (independent of \\Y\\), and
\\\delta\\ is re-drawn by the cell of the partition they induce where
\\y^\*\\ falls: for an upper bound \\c\\, \\y^\* \le c\\ gives \\\delta
= 1\\ on \\\[\epsilon, c\]\\, otherwise \\\delta = 2\\ on \\\[c, 1 -
\epsilon\]\\ (an interval \\\[l, u\]\\ gives three cells). When the
original design had more thresholds than a row records (e.g. Mode 4
intervals cut from a finer instrument) the bootstrap is conservative,
and for Mode 2 rows with a forced \\\delta\\ (a threshold built from the
score itself, which is informative) it is only an approximation. The
response must be a variable (not an expression such as `I(y / 10)`).

Replicates that fail (refit error, non-convergence, non-finite
estimates) are discarded and counted: attributes `"n_failed"` and
`"fail_rate"`, also printed. Intervals are computed from the bootstrap
distribution of each parameter, with the method controlled by `ci_type`:
`"percentile"` (default) uses the raw empirical quantiles; `"basic"`
uses reflected empirical quantiles; `"normal"` uses a normal
approximation from the bootstrap standard error (no quantiles); `"bca"`
uses bias-corrected-and-accelerated adjusted quantiles. With parametric
resampling and a nonparametric (leave-one-out) jackknife acceleration,
`"bca"` is an approximation; a warning says so once per session, and
failed jackknife refits are reported in `"n_jack_failed"`.

## Methods (by generic)

- `print(brs_bootstrap)`: Print method for bootstrap results

## Cost of `ci_type = "bca"`

The bias-corrected and accelerated interval needs an acceleration
constant, which is obtained here by a leave-one-out jackknife. That
requires `n` additional model fits, one per observation, on top of the
`R` bootstrap replicates: the total is `R + n` fits rather than `R`. The
cost is therefore driven by the sample size, not by `R`, and grows
quickly – for `n = 1000` the jackknife alone dominates the run time by
an order of magnitude. The other three interval types need only the `R`
replicates. Prefer `"percentile"` or `"basic"` for exploratory work on
large samples, and reserve `"bca"` for a final result.

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
#>   Failed replicates: 0 (0.0% of attempts)
#> 
#>     parameter   estimate   se_boot   ci_lower   ci_upper mcse_lower mcse_upper
#> 1 (Intercept)  0.2551000 0.7772708 -0.8231021 1.93714711  0.1651500 0.15824662
#> 2          x1 -0.2202060 0.4898420 -1.3384456 0.50989157  0.1751015 0.08461005
#> 3       (phi) -0.3929144 0.3259027 -1.2359624 0.05417615  0.1536645 0.12221656
#>   wald_lower wald_upper level
#> 1 -1.4390767  1.9492767  0.95
#> 2 -1.2809286  0.8405165  0.95
#> 3 -0.9343775  0.1485488  0.95
# }
```
