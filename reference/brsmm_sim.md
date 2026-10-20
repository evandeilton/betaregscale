# Simulate clustered scores from a mixed beta interval model

Simulates scores on \\0, 1, \ldots, K\\ (\\K =\\ `ncuts`) from the model
of
[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md):
Gaussian random effects in the linear predictor of the first parameter,
a latent beta value per row, and the coarsening of
[`brs_sim`](https://evandeilton.github.io/betaregscale/reference/brs_sim.md).
The output can be passed to
[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md)
as it is.

## Usage

``` r
brsmm_sim(
  formula,
  random = ~1 | id,
  data,
  beta,
  phi = 1/5,
  zeta = NULL,
  re_sd,
  re_cor = NULL,
  link = NULL,
  link_phi = NULL,
  ncuts = 100L,
  lim = 0.5,
  repar = 2L,
  interval = c("mid", "right", "left"),
  cutpoints = NULL
)
```

## Arguments

- formula:

  Model formula with one (mean) or two parts (mean `|` dispersion); a
  left-hand side is ignored.

- random:

  Random-effects formula `~ terms | group`, as in
  [`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md).

- data:

  Data frame with the covariates and the grouping variable.

- beta:

  Fixed-effect coefficients of the first parameter, on its link scale.

- phi:

  Second parameter on its link scale (one-part formulas).

- zeta:

  Coefficients of the second parameter, on its link scale (two-part
  formulas).

- re_sd:

  Standard deviations of the random effects, one per column of the
  random-effects design (e.g. intercept, then slope); positive.

- re_cor:

  Correlation matrix of the random effects; `NULL` (default) for
  independent effects. With two random effects a single correlation is
  accepted.

- link, link_phi:

  Links; `NULL` (default) selects those implied by `repar` (see
  [`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md)).

- ncuts:

  Integer \\K\\: the maximum score (default 100).

- lim:

  Half-width of the cell under `"mid"` (default 0.5).

- repar:

  Parameterization (0, 1 or 2); see
  [`brs_repar`](https://evandeilton.github.io/betaregscale/reference/brs_repar.md).

- interval:

  `"mid"` (default), `"right"` or `"left"`.

- cutpoints:

  `NULL` (default) or `c(c_lo, c_hi)`, boundary cutpoints on \\(0, 1)\\,
  as in
  [`brs_sim`](https://evandeilton.github.io/betaregscale/reference/brs_sim.md).

## Value

`data` with the columns `left`, `right`, `yt`, `y` (score) and `delta`
in front (existing columns with these names are replaced), and the
attributes `"is_prepared"`, `"ncuts"`, `"lim"` and `"interval"` that
[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md)
reuses. Attribute `"ranef"`: the simulated \\b_g\\, one row per group
level. Attribute `"truth"`: a list with `beta`, `phi` or `zeta`,
`re_sd`, `re_cor`, `D`, `L`, `cutpoints` and `par`, the true parameters
in the order and on the scale of
[`coef()`](https://rdrr.io/r/stats/coef.html) of a `brsmm` fit.

## Details

For group \\g\\ of the grouping variable of `random`, \\b_g \sim N(0,
D)\\ with \\D = \mathrm{diag}(s)\\ R\\ \mathrm{diag}(s)\\, \\s =\\
`re_sd` and \\R =\\ `re_cor` (the standard deviations and correlations
reported by
[`summary.brsmm`](https://evandeilton.github.io/betaregscale/reference/summary.brsmm.md)).
Row \\i\\ of group \\g\\ has \$\$g_1(\mu\_{gi}) = x\_{gi}^\top \beta +
x\_{r,gi}^\top b_g, \qquad g_2(\phi\_{gi}) = z\_{gi}^\top \zeta,\$\$ a
latent \\Y\_{gi} \sim \mathrm{Beta}\\ with the shapes of
[`brs_repar`](https://evandeilton.github.io/betaregscale/reference/brs_repar.md),
and the score of
[`brs_sim`](https://evandeilton.github.io/betaregscale/reference/brs_sim.md):
\\\mathrm{round}(K Y)\\ under `"mid"`, \\\lfloor (K + 1) Y \rfloor\\
under `"right"`/`"left"`, with the boundary cutpoints `c(c_lo, c_hi)`
when given. Score 0 is left-censored, \\K\\ right-censored and the rest
interval-censored.

The random effects are drawn first (a \\G \times q\\ matrix of standard
normals, \\G\\ groups in the order of `levels(factor(group))`, times
\\L^\top\\, with \\D = LL^\top\\ the lower Cholesky factor), then the
\\n\\ beta values; [`set.seed()`](https://rdrr.io/r/base/Random.html)
before the call makes the output reproducible, and
[`simulate.brsmm`](https://evandeilton.github.io/betaregscale/reference/simulate.brs.md)
uses the same sequence.

## See also

[`brs_sim`](https://evandeilton.github.io/betaregscale/reference/brs_sim.md),
[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md),
[`simulate.brsmm`](https://evandeilton.github.io/betaregscale/reference/simulate.brs.md)

## Examples

``` r
# 100 patients with 5 visits; random intercept and time slope
d <- expand.grid(visit = 1:5, id = 1:100)
d$tc <- (d$visit - 3) / 2
set.seed(7)
s <- brsmm_sim(~ tc, random = ~ 1 + tc | id, data = d,
               beta = c(-0.5, 0.4), phi = qlogis(0.2),
               re_sd = c(0.8, 0.3), re_cor = 0.3, ncuts = 10)
head(s)
#>      left   right      yt  y delta visit id   tc
#> 1 0.55000 0.65000 0.60000  6     3     1  1 -1.0
#> 2 0.45000 0.55000 0.50000  5     3     2  1 -0.5
#> 3 0.95000 0.99999 0.99999 10     2     3  1  0.0
#> 4 0.85000 0.95000 0.90000  9     3     4  1  0.5
#> 5 0.85000 0.95000 0.90000  9     3     5  1  1.0
#> 6 0.00001 0.05000 0.00001  0     1     1  2 -1.0
head(attr(s, "ranef"))
#>   (Intercept)            tc
#> 1   1.8297977  0.3543972813
#> 2  -0.9574173  0.0604336957
#> 3  -0.5554340 -0.0851900010
#> 4  -0.3298344 -0.3731870682
#> 5  -0.7765387  0.0009900384
#> 6  -0.7578240 -0.5442559799
attr(s, "truth")$par
#>                    (Intercept)                             tc 
#>                     -0.5000000                      0.4000000 
#>              (phi)_(Intercept) (re_chol_logsd)_(Intercept)|id 
#>                     -1.3862944                     -0.2231436 
#>    (re_chol)_tc:(Intercept)|id          (re_chol_logsd)_tc|id 
#>                      0.0900000                     -1.2511281 

# Random intercept only, boundary cutpoints (heaping at 0 and 10)
set.seed(8)
s1 <- brsmm_sim(~ tc, random = ~ 1 | id, data = d, beta = c(-0.5, 0.4),
                phi = qlogis(0.2), re_sd = 0.8, ncuts = 10,
                cutpoints = c(0.1, 0.9))
table(s1$y)
#> 
#>   0   1   2   3   4   5   6   7   8   9  10 
#> 105  38  66  52  50  50  35  44  33  17  10 
# \donttest{
m <- brsmm(y ~ tc, random = ~ 1 | id, data = s1, cutpoints = "boundary")
round(cbind(true = attr(s1, "truth")$par, est = coef(m)), 2)
#>                                 true   est
#> (Intercept)                    -0.50 -0.63
#> tc                              0.40  0.45
#> (phi)_(Intercept)              -1.39 -1.33
#> (re_chol_logsd)_(Intercept)|id -0.22 -0.18
#> (cut_lo)                        0.69  0.61
#> (cut_hi)                       -0.69  0.22
# }
```
