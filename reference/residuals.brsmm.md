# Residuals from a brsmm model

Residuals of a `"brsmm"` fit. The types of
[`residuals.brs`](https://evandeilton.github.io/betaregscale/reference/residuals.brs.md)
are evaluated at the conditional modes of the random effects
([`ranef.brsmm`](https://evandeilton.github.io/betaregscale/reference/ranef.brsmm.md));
for `"rqr"` and `"cell_deviance"`, `re` selects how the random effects
enter.

## Usage

``` r
# S3 method for class 'brsmm'
residuals(
  object,
  type = c("response", "pearson", "deviance", "rqr", "weighted", "sweighted",
    "cell_deviance"),
  re = c("conditional", "marginal", "posterior"),
  n_draws = 200L,
  ...
)
```

## Arguments

- object:

  A fitted `"brsmm"` object.

- type:

  Character: `"response"` (default), `"pearson"`, `"deviance"`, `"rqr"`,
  `"weighted"`, `"sweighted"` or `"cell_deviance"` (see
  [`residuals.brs`](https://evandeilton.github.io/betaregscale/reference/residuals.brs.md)).

- re:

  `"conditional"` (default), `"marginal"` or `"posterior"`; for
  `type = "rqr"` and `"cell_deviance"` only (see Details).

- n_draws:

  Number of importance-sampling proposals per group for
  `re = "posterior"` (default 200; 1 keeps the Gaussian draw).

- ...:

  Currently ignored.

## Value

Numeric vector, one value per observation.

## Details

`re` applies to `type = "rqr"` and `"cell_deviance"`:

- `"conditional"`:

  (default) at the conditional modes \\\hat b_g\\. The modes are shrunk
  toward zero, so the `"rqr"` residuals are not exactly \\N(0, 1)\\
  under the model; the shrinkage is strongest for small groups.

- `"marginal"`:

  from the marginal (population-averaged) distribution of each
  observation, with the random effects integrated out: \\v_i \sim
  U(F_m(l_i), F_m(u_i))\\, \\F_m(y) = \int F(y \mid b)\\ \varphi(b; 0,
  \hat D)\\ db\\, computed on the quadrature rule of
  [`brs_marginal`](https://evandeilton.github.io/betaregscale/reference/brs_marginal.md)
  (for score data \\F_m(l_i)\\ and \\F_m(u_i)\\ are the marginal cdf of
  the score at \\s_i - 1\\ and \\s_i\\). The rule integrates \\z =
  u/\mathrm{SD}(u)\\ over \\\[-8.5, 8.5\]\\, so marginal probabilities
  below about \\10^{-15}\\ are approximate (this matters for
  `"cell_deviance"` in extreme tails only). Each `"rqr"` residual is
  \\N(0, 1)\\ under the model; residuals of the same group are
  correlated, since they share the group's random effects.

- `"posterior"`:

  at one draw of each \\b_g\\ from its posterior given the group's data,
  by sampling-importance-resampling: `n_draws` proposals from \\N(\hat
  b_g, H_g^{-1})\\ (\\H_g\\, the curvature of the log posterior at the
  mode, by central differences; the prior \\N(0, \hat D)\\ shape when it
  is not positive definite) weighted by posterior over proposal density,
  and one of them kept. At the true parameters, \\(b_g, y_g)\\ with
  \\b_g\\ drawn from its posterior has the joint law of the model, so
  the `"rqr"` residuals are independent \\N(0, 1)\\ (up to the
  resampling approximation, of order `1/n_draws`) without the shrinkage
  of the modes; at the estimates this holds approximately. Unlike the
  marginal residuals, the residuals of a group are then nearly
  independent.

`"rqr"` residuals use \\v_i\\ uniform on \\(F(l_i), F(u_i))\\ (survival
scale above the median, to keep the precision of upper-tail cells), kept
in \\\[10^{-10}, 1 - 10^{-10}\]\\. `"cell_deviance"` residuals are
\\\mathrm{sign}\_i \sqrt{-2 \log P_i}\\ with the conditional, marginal
or posterior-draw cell probability \\P_i\\ and expected score \\\sum_s s
P(S = s) / \sum_s P(S = s)\\, as in
[`residuals.brs`](https://evandeilton.github.io/betaregscale/reference/residuals.brs.md).
The `"rqr"` residuals and the `"posterior"` draws are random: set a seed
to reproduce them.

## See also

[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md),
[`fitted.brsmm`](https://evandeilton.github.io/betaregscale/reference/fitted.brsmm.md),
[`plot.brsmm`](https://evandeilton.github.io/betaregscale/reference/plot.brsmm.md),
[`residuals.brs`](https://evandeilton.github.io/betaregscale/reference/residuals.brs.md)

## Examples

``` r
# \donttest{
dat <- data.frame(
  y = c(
    0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
    10, 40, 55, 70, 85, 25, 35, 65, 80, 15
  ),
  x1 = rep(c(1, 2), 10),
  id = factor(rep(1:4, each = 5))
)
prep <- brs_prep(dat, ncuts = 100)
#> brs_prep: n = 20 | exact = 0, left = 1, right = 1, interval = 18
fit <- brsmm(y ~ x1, random = ~ 1 | id, data = prep)
head(residuals(fit))
#> [1] -0.3855953 -0.2593261 -0.1856053  0.1906739  0.3643947  0.3286716
head(residuals(fit, type = "pearson"))
#> [1] -1.3225376 -0.9366406 -0.6366000  0.6886809  1.2498228  1.1087368
# }
# Quantile residuals without the shrinkage of the conditional modes
d <- expand.grid(visit = 1:5, id = 1:60)
set.seed(3)
s <- brsmm_sim(~ visit, random = ~ 1 | id, data = d, beta = c(-0.6, 0.1),
               phi = qlogis(0.2), re_sd = 0.7, ncuts = 10)
m <- brsmm(y ~ visit, random = ~ 1 | id, data = s)
set.seed(1)
r_m <- residuals(m, type = "rqr", re = "marginal")
r_p <- residuals(m, type = "rqr", re = "posterior")
round(c(mean_m = mean(r_m), sd_m = sd(r_m), mean_p = mean(r_p), sd_p = sd(r_p)), 2)
#> mean_m   sd_m mean_p   sd_p 
#>   0.01   1.00  -0.01   1.00 
```
