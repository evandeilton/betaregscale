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

Each refit starts from the estimate of `object` and uses the compiled
Hessian (`hessian_method = "cpp"`), so replicates are cheap.

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
# Synthetic NRS-11 scores at 6h, 12h, 24h (time is a factor)
set.seed(3)
nrs <- data.frame(time = factor(rep(c("6h", "12h", "24h"), each = 40),
                                levels = c("6h", "12h", "24h")))
shp <- brs_repar(mu = plogis(-1.3 + c(0, 0.75, 0.3)[nrs$time]), phi = 0.3)
nrs$y <- round(10 * rbeta(nrow(nrs), shp$shape1, shp$shape2))
fit <- brs(y ~ time, data = nrs, ncuts = 10)

# Percentile intervals from 30 parametric replicates (use R >= 199 in practice)
set.seed(4)
bt <- brs_bootstrap(fit, R = 30)
bt
#> Bootstrap confidence intervals
#>   Level: 0.95 | CI: percentile | Successful replicates: 30 / 30 | Attempts: 30 
#>   Failed replicates: 0 (0.0% of attempts)
#> 
#>     parameter   estimate   se_boot    ci_lower   ci_upper mcse_lower mcse_upper
#> 1 (Intercept) -1.2953746 0.2011046 -1.71296642 -0.9501425 0.07569263 0.04900518
#> 2     time12h  0.7657720 0.2290701  0.29264721  1.1319112 0.05329737 0.09006074
#> 3     time24h  0.3854151 0.2849472 -0.06794233  0.8036531 0.04243362 0.04330853
#> 4       (phi) -0.8166662 0.1422067 -1.12895854 -0.6751089 0.03337507 0.01822655
#>   wald_lower wald_upper level
#> 1 -1.6863745 -0.9043747  0.95
#> 2  0.2513386  1.2802054  0.95
#> 3 -0.1363913  0.9072215  0.95
#> 4 -1.0759353 -0.5573970  0.95
# Bootstrap next to Wald limits, and the bootstrap/Wald SE ratio
cols <- c("parameter", "ci_lower", "ci_upper", "wald_lower", "wald_upper")
as.data.frame(bt)[, cols]
#>     parameter    ci_lower   ci_upper wald_lower wald_upper
#> 1 (Intercept) -1.71296642 -0.9501425 -1.6863745 -0.9043747
#> 2     time12h  0.29264721  1.1319112  0.2513386  1.2802054
#> 3     time24h -0.06794233  0.8036531 -0.1363913  0.9072215
#> 4       (phi) -1.12895854 -0.6751089 -1.0759353 -0.5573970
round(bt$se_boot / sqrt(diag(vcov(fit))), 2)
#> (Intercept)     time12h     time24h       (phi) 
#>        1.01        0.87        1.07        1.08 
```
