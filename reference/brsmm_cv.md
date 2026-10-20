# Group-level cross-validation for brsmm models

K-fold cross-validation of a
[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md)
model with folds formed by group: all observations of a group (cluster)
are in the same fold. Each held-out group is scored by its marginal
log-likelihood at the estimates of the fit to the other folds.

## Usage

``` r
brsmm_cv(
  formula,
  random = ~1 | id,
  data,
  k = 10L,
  repeats = 1L,
  folds = NULL,
  seed = NULL,
  score_points = NULL,
  ...
)
```

## Arguments

- formula, random:

  Model and random-effects formulas passed to
  [`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md).

- data:

  Data frame.

- k:

  Number of folds (default 10).

- repeats:

  Number of repeated k-fold runs.

- folds:

  Optional fold assignment: a vector named by the group levels (read by
  name when its unique names are exactly the levels), or else a vector
  with one fold label per row of `data`, constant within each group.
  Overrides `k`, `repeats` and `seed`.

- seed:

  Optional integer seed for the fold assignment.

- score_points:

  Nodes per dimension of the adaptive Gauss-Hermite rule that scores
  held-out groups; `NULL` for the default (see Details).

- ...:

  Further arguments for
  [`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md)
  (e.g. `ncuts`, `int_method`, `n_points`, `cutpoints`).

## Value

A data frame of class `"brsmm_cv"` with one row per fold: `repeat`,
`fold`, `n_groups_train`, `n_groups_test`, `n_train`, `n_test`,
`log_score` (sum over the held-out groups), `converged`, `error` and
`warning` (messages of the fold fit). Attributes: `log_score` (total per
repeat; `NA` if a fold failed), `group_scores` (held-out contribution of
each group, one vector per repeat), `folds` (fold of each group, one
vector per repeat, named by group level: integer for drawn folds, the
type of the supplied `folds` otherwise, character for a factor) and
`score_points`.

## Details

The held-out log score of a fold is \\\sum_g \log \int \prod_i P(S =
s\_{gi} \mid b)\\ \varphi(b; 0, \hat D)\\ db\\ over its groups, at the
training-fold estimates
([`brs_loglik_contrib`](https://evandeilton.github.io/betaregscale/reference/brs_loglik_contrib.md)
with `newdata`). The random effects of held-out groups are integrated
out, so no random-effect estimates are needed. The integral is computed
by adaptive Gauss-Hermite quadrature with `score_points` nodes per
random-effect dimension, whatever `int_method` the fits use; the default
is 31 for one or two random effects and otherwise the largest number
whose product grid has at most 20000 nodes (27 for three, 11 for four, 7
for five; at least 2). The grid may not exceed 500000 nodes. For score
data the log score is the log probability of the observed scores (no
Jacobian term), so it is on the same scale as the held-out
log-likelihood of any model that gives a probability mass function on
\\0, \ldots, K\\, evaluated on the same folds (attribute `"folds"`).

Rows with missing values in the model variables (response, covariates,
random-effect terms or group) are dropped first; a group without
complete rows is not used, and `n_train` / `n_test` count complete rows.

Folds: `k` folds of nearly equal numbers of groups, drawn at random for
each repeat. With `seed`, the folds are drawn after `set.seed(seed)` and
the RNG state is restored afterwards; without it, the current RNG state
is used. `folds` gives the assignment directly (one repeat). Held-out
rows use the cells of the training fit (its estimated cutpoints under
`cutpoints = "boundary"`). A fold whose fit or scoring fails (for
example, a factor level or a border score missing from the training
folds) has `log_score = NA` and the message in `error`.

## See also

[`brs_cv`](https://evandeilton.github.io/betaregscale/reference/brs_cv.md),
[`brs_loglik_contrib`](https://evandeilton.github.io/betaregscale/reference/brs_loglik_contrib.md),
[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md)

## Examples

``` r
# \donttest{
set.seed(21)
nrs <- expand.grid(id = 1:30, time = c("6h", "12h", "24h"))
b0 <- rnorm(30, sd = 0.8)
eta <- -1.3 + c(0, 0.75, 0.3)[nrs$time] + b0[nrs$id]
shp <- brs_repar(mu = plogis(eta), phi = 0.2, repar = 2)
nrs$y <- round(10 * rbeta(nrow(nrs), shp$shape1, shp$shape2))

cv <- brsmm_cv(y ~ time, random = ~ 1 | id, data = nrs, k = 5, seed = 1,
               ncuts = 10)
cv
#>   repeat fold n_groups_train n_groups_test n_train n_test log_score converged
#> 1      1    1             24             6      72     18 -37.84983      TRUE
#> 2      1    2             24             6      72     18 -40.70206      TRUE
#> 3      1    3             24             6      72     18 -39.15762      TRUE
#> 4      1    4             24             6      72     18 -39.61621      TRUE
#> 5      1    5             24             6      72     18 -41.23546      TRUE
#>   error warning
#> 1  <NA>    <NA>
#> 2  <NA>    <NA>
#> 3  <NA>    <NA>
#> 4  <NA>    <NA>
#> 5  <NA>    <NA>
attr(cv, "log_score")
#> [1] -198.5612
# }
```
