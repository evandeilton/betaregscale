# Log-likelihood contributions per observation or per group

The terms whose sum is the log-likelihood: one per observation for a
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md)
fit, and one per group (cluster) for a
[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md)
fit, the marginal log-likelihood of the group's observations with its
random effects integrated out. With `newdata`, the contributions of new
observations at the fitted parameters (held-out scoring).

## Usage

``` r
brs_loglik_contrib(object, newdata = NULL, int_method = NULL, n_points = NULL)
```

## Arguments

- object:

  A fitted `"brs"` or `"brsmm"` object.

- newdata:

  Optional data frame of new observations (see Details).

- int_method, n_points:

  For `brsmm` fits: integration rule (`"laplace"`, `"aghq"` or `"qmc"`)
  and number of nodes (nodes per dimension for `"aghq"`, points for
  `"qmc"`); `NULL` (default) uses those of the fit. Ignored for `brs`
  fits.

## Value

Named numeric vector: by row name of the data (`brs`) or by group level
(`brsmm`, with attribute `"n"`, the number of observations of each
group).

## Details

For a `brs` fit the contribution of a score row is \\\log P(S = s_i)\\,
the log probability of its cell (the log density for exact observations,
\\\delta = 0\\), computed as in the compiled likelihood. For a `brsmm`
fit the contribution of group \\g\\ is \$\$\log \int \prod_i P(S =
s\_{gi} \mid b)\\ \varphi(b; 0, D)\\ db,\$\$ computed by the fit's
`int_method` and `n_points` (or `qmc_points`) unless `int_method` /
`n_points` are given. With the defaults and `newdata = NULL`, the
contributions add up to `logLik(object)` (under `"laplace"`, up to the
tolerance of the inner mode search, about \\10^{-8}\\). A group without
a positive-definite random-effect mode, which adds a penalty value to
the fitted log-likelihood, gets `NA` with a warning.

With `newdata`, rows need the response, the covariates and (for `brsmm`)
the grouping variable; rows with missing values in them are dropped. The
groups need not be groups of the fit: their random effects are
integrated out, so no random-effect estimate is used. Scores are
coarsened with the cells of the fit (its estimated cutpoints under
`cutpoints = "boundary"`). For score data the contributions are log
probabilities of the observed scores (no Jacobian term), on the same
scale as the log-likelihood of any model that gives a probability mass
function on \\0, \ldots, K\\.

## See also

[`logLik.brsmm`](https://evandeilton.github.io/betaregscale/reference/logLik.brsmm.md),
[`brsmm_cv`](https://evandeilton.github.io/betaregscale/reference/brsmm_cv.md),
[`brs_marginal`](https://evandeilton.github.io/betaregscale/reference/brs_marginal.md)

## Examples

``` r
set.seed(21)
nrs <- expand.grid(id = 1:30, time = c("6h", "12h", "24h"))
b0 <- rnorm(30, sd = 0.8)
eta <- -1.3 + c(0, 0.75, 0.3)[nrs$time] + b0[nrs$id]
shp <- brs_repar(mu = plogis(eta), phi = 0.2, repar = 2)
nrs$y <- round(10 * rbeta(nrow(nrs), shp$shape1, shp$shape2))
m <- brsmm(y ~ time, random = ~ 1 | id, data = nrs, ncuts = 10)

ll_g <- brs_loglik_contrib(m)
head(ll_g)
#>         1         2         3         4         5         6 
#> -6.019919 -6.321018 -7.271889 -5.779475 -8.037281 -8.350358 
c(sum = sum(ll_g), logLik = as.numeric(logLik(m)))
#>       sum    logLik 
#> -192.9768 -192.9768 
# Same groups, adaptive Gauss-Hermite with 31 nodes
head(brs_loglik_contrib(m, int_method = "aghq", n_points = 31))
#>         1         2         3         4         5         6 
#> -6.002478 -6.304053 -7.253131 -5.763656 -8.019587 -8.333486 
```
