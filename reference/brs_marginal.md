# Population-averaged score probabilities, expected scores and exceedance

Score distribution of a
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md) or
[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md)
fit at given covariate values, with standard errors and intervals: the
probabilities \\P(S = s)\\ of the scores \\0, \ldots, K\\, the expected
score or the exceedance probabilities \\P(S \ge k)\\, and differences of
these between pairs of rows. For `brsmm` fits the random effects are
integrated out (marginal, population-averaged quantities) or set to zero
(`conditional = TRUE`).

## Usage

``` r
brs_marginal(
  object,
  newdata = NULL,
  type = c("prob", "expected_score", "exceed"),
  k = NULL,
  contrast = NULL,
  conditional = FALSE,
  level = 0.95,
  se_method = c("delta", "simulation"),
  n_sim = 1000L,
  n_points = 8L
)
```

## Arguments

- object:

  A fitted `"brs"` or `"brsmm"` object.

- newdata:

  Data frame of covariate values, one prediction per row; `NULL`
  (default) uses the rows of the fit.

- type:

  `"prob"` (default; probabilities of the scores \\0, \ldots, K\\),
  `"expected_score"` or `"exceed"`.

- k:

  Integer thresholds in \\0, \ldots, K\\ for `type = "exceed"`; `NULL`
  (default) uses \\1, \ldots, K\\.

- contrast:

  Optional pairs of rows of `newdata` to compare: a two-column matrix
  with one pair per row, a list of length-2 vectors, or a single
  length-2 vector. The pair `c(i, j)` gives the quantity of row `i`
  minus that of row `j`.

- conditional:

  Logical; `TRUE` sets the random effects to zero instead of integrating
  them out (no effect for `brs` fits).

- level:

  Confidence level of the intervals (default 0.95).

- se_method:

  `"delta"` (default) or `"simulation"`.

- n_sim:

  Number of parameter draws for `se_method = "simulation"` (default
  1000).

- n_points:

  Gauss-Legendre nodes per panel of the random-effect integral (default
  8; see Details).

## Value

A data frame of class `"brs_marginal"` with the columns of `newdata`
(or, for contrasts, `contrast`, `row1` and `row2`), `row`, `score`
(`type = "prob"`) or `k` (`type = "exceed"`), `estimate`, `std.error`,
`lower` and `upper`. Attributes: `type`, `conditional`, `level`,
`se_method`, `n_nodes` (total quadrature nodes; `NA` without random
effects) and, for the delta method, `jacobian` (quantities by
parameters).

## Details

For a row with fixed-effect linear predictor \\\eta\\ and random-effect
design \\x_r\\, the marginal score probabilities are \$\$P(S = s) = \int
P(S = s \mid \eta + x_r^\top b)\\ \varphi(b; 0, D)\\ db = E\\P(S = s
\mid \eta + u)\\, \qquad u \sim N(0, x_r^\top D x_r),\$\$ because \\b\\
enters the model only through \\u = x_r^\top b\\ (the dispersion
submodel has no random effects), so the integral is one-dimensional
whatever the number of random effects. It is computed on \\z = u /
\mathrm{SD}(u) \in \[-8.5, 8.5\]\\ by a composite Gauss-Legendre rule
with `n_points` nodes per panel and panels of width at most 1 on the
scale of \\z\\ and at most 0.5 on the scale of \\u\\ (the conditional
probabilities can change quickly with \\u\\ when the precision is high
or the variance large, where Gauss-Hermite rules converge slowly). A
warning is given when halving the panel width changes a value by more
than \\10^{-6}\\.

\\P(S = s \mid \cdot)\\ is the probability of the cell of \\s\\ under
the latent beta law, with the cells of the fit's `interval` map or its
estimated boundary cutpoints (as in
[`brs_predict_scoreprob`](https://evandeilton.github.io/betaregscale/reference/brs_predict_scoreprob.md),
but without a probability floor), so each row sums to 1 whenever the
cells partition (0, 1) (`lim = 0.5`, or `interval = "right"` or
`"left"`). `"expected_score"` is \\\sum_s s\\ P(S = s)\\ on the \\0,
\ldots, K\\ scale and `"exceed"` is \\P(S \ge k) = \sum\_{s \ge k} P(S =
s)\\. When the cells do not partition (0, 1) (`interval = "mid"` with
`lim < 0.5`), the probabilities sum to less than 1 and these two types
are not moments or tail probabilities of a distribution; they are
computed with a warning.

Under a nonlinear link the marginal quantities differ from the
conditional ones at \\b = 0\\: the latter describe a subject whose
random effects are zero, not the population average, and the gap grows
with the random-effect variance. For a `brs` fit both coincide.

Uncertainty: with `se_method = "delta"` (default), delta-method standard
errors from a Richardson central-difference Jacobian of the quantities
in all parameters of `coef(object)` (regression, precision,
random-effect Cholesky and boundary-cutpoint parameters) and
`vcov(object)`. The Wald intervals are computed on the logit scale for
probabilities and for \\E\[S\]/K\\ and transformed back (so they stay in
range), and on the natural scale for contrasts. With
`se_method = "simulation"`, `n_sim` parameter vectors are drawn from
\\N(\hat\theta, V)\\; the standard error is the standard deviation of
the quantity over the draws and the interval its percentile interval
(the draws use the RNG: call
[`set.seed()`](https://rdrr.io/r/base/Random.html) first for
reproducible values). The estimate is always the plug-in value at
\\\hat\theta\\.

Rows of `newdata` with missing covariates give `NA`; the response and
the grouping variable are not needed.

## See also

[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md),
[`brs_predict_scoreprob`](https://evandeilton.github.io/betaregscale/reference/brs_predict_scoreprob.md),
[`predict.brsmm`](https://evandeilton.github.io/betaregscale/reference/predict.brsmm.md),
[`brs_loglik_contrib`](https://evandeilton.github.io/betaregscale/reference/brs_loglik_contrib.md)

## Examples

``` r
# Synthetic NRS-11 scores of 40 patients at three times. Simulated data.
set.seed(21)
nrs <- expand.grid(id = 1:40, time = c("6h", "12h", "24h"))
b0 <- rnorm(40, sd = 0.8)
eta <- -1.3 + c(0, 0.75, 0.3)[nrs$time] + b0[nrs$id]
shp <- brs_repar(mu = plogis(eta), phi = 0.2, repar = 2)
nrs$y <- round(10 * rbeta(nrow(nrs), shp$shape1, shp$shape2))
m <- brsmm(y ~ time, random = ~ 1 | id, data = nrs, ncuts = 10)

nd <- data.frame(time = c("6h", "12h", "24h"))
# Population-averaged mean score and P(S >= 4) at each time
brs_marginal(m, nd, type = "expected_score")
#>   time row estimate std.error    lower    upper
#> 1   6h   1 2.765696 0.3641913 2.111009 3.532505
#> 2  12h   2 3.851484 0.4053020 3.093409 4.669728
#> 3  24h   3 3.270772 0.3888986 2.558317 4.073042
brs_marginal(m, nd, type = "exceed", k = 4)
#>   time row k  estimate  std.error     lower     upper
#> 1   6h   1 4 0.3343305 0.05886989 0.2302123 0.4575491
#> 2  12h   2 4 0.5081511 0.06331149 0.3860663 0.6292714
#> 3  24h   3 4 0.4157595 0.06246979 0.3006367 0.5408728
# Same at b = 0 (a typical patient): not the population average
brs_marginal(m, nd, type = "expected_score", conditional = TRUE)
#>   time row estimate std.error    lower    upper
#> 1   6h   1 2.505482 0.3874073 1.824162 3.337397
#> 2  12h   2 3.693577 0.4560519 2.852210 4.622625
#> 3  24h   3 3.050064 0.4262106 2.283506 3.942478
# Change from 6h: 12h - 6h and 24h - 6h
brs_marginal(m, nd, type = "expected_score",
             contrast = rbind(c(2, 1), c(3, 1)))
#>   contrast row1 row2  estimate std.error      lower    upper
#> 1    2 - 1    2    1 1.0857885 0.4060150  0.2900137 1.881563
#> 2    3 - 1    3    1 0.5050766 0.3970665 -0.2731594 1.283313
```
