# Simulate responses from a fitted beta interval model

Draws new responses from a
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md) or
[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md)
fit at the estimated parameters, coarsened with the cells of the fit, as
[`brs_sim`](https://evandeilton.github.io/betaregscale/reference/brs_sim.md)
and
[`brsmm_sim`](https://evandeilton.github.io/betaregscale/reference/brsmm_sim.md)
do.

## Usage

``` r
# S3 method for class 'brs'
simulate(object, nsim = 1, seed = NULL, ...)

# S3 method for class 'brsmm'
simulate(object, nsim = 1, seed = NULL, re = c("marginal", "conditional"), ...)
```

## Arguments

- object:

  A fitted `"brs"` or `"brsmm"` object.

- nsim:

  Number of simulated response vectors (default 1).

- seed:

  `NULL` or an integer seed, as in
  [`simulate`](https://rdrr.io/r/stats/simulate.html): with a seed the
  caller's RNG state is restored afterwards.

- ...:

  Currently ignored.

- re:

  For `brsmm` fits: `"marginal"` (default; new random effects in each
  simulation) or `"conditional"` (random effects at the conditional
  modes).

## Value

A data frame with `nsim` columns `sim_1`, ..., one row per observation
of the fit, and the attribute `"seed"` (as
[`simulate`](https://rdrr.io/r/stats/simulate.html)). For `brsmm` fits
with `re = "marginal"`, the attribute `"ranef"` holds the random effects
of each simulation (a list of matrices, one row per group).

## Details

Each simulation draws one latent beta value per row with the fitted
shapes and returns its score on \\0, \ldots, K\\ (cells of the fit's
`interval` map, or its estimated boundary cutpoints); rows fitted as
exact values (\\\delta = 0\\) get the latent value itself.

For a `brsmm` fit, `re = "marginal"` (default) draws new random effects
\\b_g \sim N(0, \hat D)\\ for each simulation (new groups with the same
design: the marginal distribution of the data), in the order of
[`brsmm_sim`](https://evandeilton.github.io/betaregscale/reference/brsmm_sim.md),
so `simulate(fit, seed = s)` and `set.seed(s); brsmm_sim(...)` with the
estimates as true values give the same scores. `re = "conditional"`
keeps the random effects at their conditional modes
([`ranef.brsmm`](https://evandeilton.github.io/betaregscale/reference/ranef.brsmm.md))
and draws only the latent values: the same groups again.

## See also

[`brs_sim`](https://evandeilton.github.io/betaregscale/reference/brs_sim.md),
[`brsmm_sim`](https://evandeilton.github.io/betaregscale/reference/brsmm_sim.md)

## Examples

``` r
set.seed(2023)
d <- data.frame(x = runif(150))
s <- brs_sim(~ x, data = d, beta = c(-0.5, 1), phi = qlogis(0.3), ncuts = 10)
fit <- brs(y ~ x, data = s)
sims <- simulate(fit, nsim = 200, seed = 1)
# Observed and simulated frequencies of the scores
rbind(observed = tabulate(s$y + 1, 11) / nrow(s),
      simulated = tabulate(unlist(sims) + 1, 11) / (200 * nrow(s)))
#>             [,1]       [,2]      [,3]      [,4]   [,5]      [,6]      [,7]
#> observed  0.0400 0.07333333 0.1066667 0.1266667 0.1200 0.1266667 0.1200000
#> simulated 0.0418 0.09500000 0.1001333 0.1044667 0.1065 0.1078333 0.1032333
#>                 [,8]       [,9]      [,10]      [,11]
#> observed  0.06666667 0.08666667 0.07333333 0.06000000
#> simulated 0.10476667 0.10166667 0.09416667 0.04043333
```
