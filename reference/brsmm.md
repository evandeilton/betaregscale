# Fit a mixed-effects beta interval regression model

Beta interval regression
([`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md))
with Gaussian random effects in the linear predictor of the first
parameter (the mean, or the shape \\p\\ under `repar = 0`), fitted by
marginal maximum likelihood. `random = ~ 1 | id` gives a random
intercept per group, `~ 1 + x | id` a random intercept and slope with a
free correlation.

## Usage

``` r
brsmm(
  formula,
  random = ~1 | id,
  data,
  link = NULL,
  link_phi = NULL,
  repar = 2L,
  ncuts = NULL,
  lim = NULL,
  int_method = c("laplace", "aghq", "qmc"),
  n_points = 11L,
  qmc_points = 1024L,
  start = NULL,
  method = c("BFGS", "L-BFGS-B"),
  hessian_method = c("cpp", "numDeriv", "optim"),
  control = list(maxit = 2000L),
  interval = NULL
)
```

## Arguments

- formula:

  Model formula: `y ~ x1 + x2` or `y ~ x1 + x2 | z1 + z2` (see
  [`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md)).

- random:

  Random-effects formula `~ terms | group`, e.g. `~ 1 | id` or
  `~ 1 + x | id`.

- data:

  Data frame (raw scores, or the output of
  [`brs_prep`](https://evandeilton.github.io/betaregscale/reference/brs_prep.md)).

- link, link_phi:

  Links for the first and second parameter; `NULL` (default) selects
  those implied by `repar` (see
  [`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md)).

- repar:

  Parameterisation (0, 1 or 2); see
  [`brs_repar`](https://evandeilton.github.io/betaregscale/reference/brs_repar.md).

- ncuts:

  Integer \\K\\: the maximum score (scale \\0, \ldots, K\\). `NULL`
  (default) uses the value stored by
  [`brs_prep`](https://evandeilton.github.io/betaregscale/reference/brs_prep.md),
  or 100; an explicit different value is ignored with a warning.

- lim:

  Half-width of the cell in \\(0, 0.5\]\\ (`interval = "mid"` only);
  `NULL` uses the stored value, or 0.5.

- int_method:

  `"laplace"` (default), `"aghq"` or `"qmc"`; see Details.

- n_points:

  Nodes per dimension for `"aghq"` (default 11).

- qmc_points:

  Halton points for `"qmc"` (default 1024).

- start:

  Optional starting vector: fixed effects, precision coefficients, then
  the packed Cholesky parameters.

- method:

  `"BFGS"` (default) or `"L-BFGS-B"`.

- hessian_method:

  `"cpp"` (default; Richardson differences of the compiled gradient),
  `"numDeriv"` or `"optim"`.

- control:

  Control list for [`optim`](https://rdrr.io/r/stats/optim.html), merged
  into the default `list(maxit = 2000L)`.

- interval:

  `"mid"`, `"right"` or `"left"` (see
  [`brs_check`](https://evandeilton.github.io/betaregscale/reference/brs_check.md));
  `NULL` uses the stored value, or `"mid"`.

## Value

An object of class `"brsmm"` with the components of a `"brs"` fit
(`par`, `coefficients` with a `random` part, `value`, `hessian`,
`diagnostics`, ...) and `random` (group variable, levels, conditional
modes `mode_b`, `D`, `L` and the SDs `sd_b`), `ngroups`, `int_method`.

## Details

For group \\g\\ with observations \\i\\, the model is \$\$g_1(\mu\_{gi})
= x\_{gi}^\top \beta + x\_{r,gi}^\top b_g, \qquad g_2(\phi\_{gi}) =
z\_{gi}^\top \gamma, \qquad b_g \sim N(0, D),\$\$ with the censored
contributions of
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md).
The marginal log-likelihood is \\\sum_g \log \int \exp\\h_g(b)\\\\ db\\,
where \\h_g(b) = \sum_i \ell\_{gi}(b) + \log \varphi(b; 0, D)\\. The
integral is computed around the mode \\\hat b_g\\ of \\h_g\\, with \\H_g
= -\partial^2 h_g / \partial b\\ \partial b^\top\\ there:

- `"laplace"`:

  \\h_g(\hat b_g) + \frac{q}{2}\log(2\pi) - \frac12 \log \|H_g\|\\;
  fast, accurate when groups are not tiny.

- `"aghq"`:

  adaptive Gauss-Hermite quadrature: a product grid of `n_points` nodes
  per dimension at \\\hat b_g + \sqrt{2}\\ H_g^{-1/2} z\\, with the
  symmetric square root \\H_g^{-1/2}\\ (`n_points`\\^q\\ nodes, at most
  500000).

- `"qmc"`:

  importance sampling from \\N(\hat b_g, H_g^{-1})\\ (nodes \\\hat b_g +
  H_g^{-1/2} z\\) on `qmc_points` Halton points. It is deterministic
  and, with two or more random effects, underestimates the
  log-likelihood (about \\-0.05\\ at 1024 points in the package's
  two-effect checks); prefer `"aghq"` for up to three random effects.

The inner mode is found by a Levenberg–Marquardt Newton method,
warm-started from the modes of the previous evaluation (the cache is
cleared at the start of each fit). A group whose curvature is not
positive definite at its mode adds the penalty value \\-10^6\\ instead
of a silently regularised term, and `brsmm()` warns.

\\D = LL^\top\\ is parameterised by its lower Cholesky factor: the log
of each diagonal entry and the off-diagonal entries as they are
(`(re_chol_logsd)_` and `(re_chol)_` in
[`coef()`](https://rdrr.io/r/stats/coef.html)).
[`summary.brsmm`](https://evandeilton.github.io/betaregscale/reference/summary.brsmm.md)
reports the standard deviations and correlations with intervals;
[`brsmm_re_study`](https://evandeilton.github.io/betaregscale/reference/brsmm_re_study.md)
the intraclass correlation.

## Estimation

[`optim`](https://rdrr.io/r/stats/optim.html) maximises the marginal
log-likelihood with its compiled gradient: the derivative of the chosen
approximation by the chain rule on the linear predictor and the
implicit-function theorem at the group modes (including the movement of
the quadrature nodes), with per-observation derivatives by central
differences. The Hessian for
[`vcov()`](https://rdrr.io/r/stats/vcov.html) (`hessian_method = "cpp"`,
the default) is a Richardson central difference of that gradient.
Starting values: `start` when given; otherwise those of
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md) for
the fixed effects and \\\log\\ of the between-group SD of the cell
centres (at least 0.1) for the random effects.

## Diagnostics

The checks of
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md)
('Fit diagnostics') apply, with the compiled gradient (a central
difference with step \\10^{-3}\\ if it is not finite). In addition:

- “Variance component at the boundary”:

  A random-effect term has log SD below \\-6\\, or raises the
  log-likelihood by less than \\10^{-3}\\ over the same fit with that SD
  at zero. Its variance is essentially zero; the standard error of its
  log SD is meaningless. Test the term with
  [`anova.brsmm`](https://evandeilton.github.io/betaregscale/reference/anova.brsmm.md)
  (chi-bar-square mixture) and drop it if not needed. When the gain is
  negative the message adds that SD \\\approx 0\\ has a higher
  log-likelihood: `optim` stopped short of the maximum.

- “group(s) have no positive-definite random-effect mode”:

  Those groups contribute the penalty value; the fit is not reliable.
  Simplify the random-effects structure or check the data of those
  groups.

`fit$diagnostics` stores `re_boundary`, `re_gain` and `inner` (number of
such groups, largest gradient norm at the modes). Rank-deficient
fixed-effect or random-effect design matrices are an error.

## References

Lopes, J. E. (2023). *Modelos de regressao beta para dados de escala*.
Master's dissertation, Universidade Federal do Parana, Curitiba. URI:
https://hdl.handle.net/1884/86624.

Ferrari, S. L. P., and Cribari-Neto, F. (2004). Beta regression for
modelling rates and proportions. *Journal of Applied Statistics*,
**31**(7), 799–815.
[doi:10.1080/0266476042000214501](https://doi.org/10.1080/0266476042000214501)

Pinheiro, J. C., and Bates, D. M. (1995). Approximations to the
log-likelihood function in the nonlinear mixed-effects model. *Journal
of Computational and Graphical Statistics*, **4**(1), 12–35.
[doi:10.1080/10618600.1995.10474663](https://doi.org/10.1080/10618600.1995.10474663)

## See also

[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md),
[`summary.brsmm`](https://evandeilton.github.io/betaregscale/reference/summary.brsmm.md),
[`anova.brsmm`](https://evandeilton.github.io/betaregscale/reference/anova.brsmm.md),
[`brsmm_re_study`](https://evandeilton.github.io/betaregscale/reference/brsmm_re_study.md),
[`ranef.brsmm`](https://evandeilton.github.io/betaregscale/reference/ranef.brsmm.md),
[`predict.brsmm`](https://evandeilton.github.io/betaregscale/reference/predict.brsmm.md)

## Examples

``` r
# Synthetic NRS-11 scores (0-10) of 40 patients at 6h, 12h and 24h;
# intercepts and time slopes vary by patient. Simulated, not real data.
set.seed(21)
nrs <- expand.grid(id = 1:40, time = c("6h", "12h", "24h"))
nrs$tc <- c(-1, 0, 1)[nrs$time]                 # centred time for the slope
b0 <- rnorm(40, sd = 0.8)
b1 <- rnorm(40, sd = 0.5)
eta <- -1.3 + c(0, 0.75, 0.3)[nrs$time] + b0[nrs$id] + b1[nrs$id] * nrs$tc
shp <- brs_repar(mu = plogis(eta), phi = 0.2, repar = 2)
nrs$y <- round(10 * rbeta(nrow(nrs), shp$shape1, shp$shape2))

# Random intercept per patient (Laplace): SD with its interval, no z-test
m1 <- brsmm(y ~ time, random = ~ 1 | id, data = nrs, ncuts = 10)
summary(m1)
#> 
#> Call:
#> brsmm(formula = y ~ time, random = ~1 | id, data = nrs, ncuts = 10)
#> 
#> Randomized Quantile Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -2.4305 -0.6884  0.0513  0.6250  2.3968 
#> 
#> Coefficients (mean model with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  -1.1663     0.2094  -5.570 2.55e-08 ***
#> time12h       0.7683     0.2472   3.108  0.00188 ** 
#> time24h       0.3025     0.2515   1.203  0.22896    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Phi coefficients (precision model with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  -1.0747     0.1715  -6.266 3.71e-10 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Random effects (SD and Corr; 95% Wald CI on the log / atanh scale; no tests, see anova()):
#>                Estimate Lower  Upper
#> SD (Intercept)   0.5607 0.315 0.9983
#> ---
#> Mixed beta interval model (Laplace)
#> Observations: 120  | Groups: 40 
#> Log-likelihood: -262.0285 on 5 Df | AIC: 534.0570 | BIC: 547.9945 
#> Pseudo R-squared: 0.0669 
#> Number of iterations: 28 (BFGS) 
#> Censoring: 96 interval | 22 left | 2 right 
#> 

# Random intercept and slope: SDs and their correlation
m2 <- brsmm(y ~ time, random = ~ 1 + tc | id, data = nrs, ncuts = 10)
summary(m2)$varcorr
#>                  term type  estimate      lower     upper se_transformed
#> 1      SD (Intercept)   sd 0.6328837  0.3894562 1.0284643      0.2477267
#> 2               SD tc   sd 0.5474364  0.2490126 1.2034998      0.4019169
#> 3 Corr tc,(Intercept) corr 0.7034330 -0.4083700 0.9748543      0.6672150
head(ranef(m2))
#>   (Intercept)         tc
#> 1   0.3540762  0.1152133
#> 2   0.1124897 -0.1982327
#> 3   1.1300445  0.6532920
#> 4  -0.2328539 -0.1960478
#> 5   0.9638750  0.8034315
#> 6   0.1412063  0.1006265

# Is the slope needed? One added random term: the p-value uses the
# mixture 1/2 chi2(1) + 1/2 chi2(2) (see the printed heading)
anova(m1, m2)
#> Likelihood-ratio comparison of brs/brsmm models
#> Rows M2: one added random effect (variance on the boundary); Pr(>Chisq) from the chi-bar-square mixture 1/2 chi2(Df - 1) + 1/2 chi2(Df).
#> 
#>            Df  logLik    AIC    BIC  Chisq Chi Df Pr(>Chisq)  
#> M1 (brsmm)  5 -262.03 534.06 547.99                           
#> M2 (brsmm)  7 -259.15 532.30 551.81 5.7547      2    0.03636 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# A known patient uses its conditional mode; a new one gets b = 0
nd <- data.frame(time = "12h", tc = 0, id = c(1, 999))
predict(m2, newdata = nd)
#> [1] 0.4840601 0.3970295
```
