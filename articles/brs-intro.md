# Introduction to betaregscale

This vignette describes the model fitted by **betaregscale**: how scores
are mapped to intervals of $`(0, 1)`$, the beta law in mean–dispersion
form, the likelihood, the choice of interval map and the estimable
boundary cutpoints. The other vignettes cover mixed models
([`vignette("brs-mm")`](https://evandeilton.github.io/betaregscale/articles/brs-mm.md)),
post-estimation tools
([`vignette("brs-analyst-tools")`](https://evandeilton.github.io/betaregscale/articles/brs-analyst-tools.md))
and longer workflows
([`vignette("brs-advanced-workflows")`](https://evandeilton.github.io/betaregscale/articles/brs-advanced-workflows.md)).

``` r

library(betaregscale)
```

## The score scale

Scores run over $`0, 1, \ldots, K`$, where $`K`$ is the **maximum
score** (argument `ncuts`). The scale has $`K + 1`$ categories: an
NRS-11 pain scale (0 to 10) has `ncuts = 10`, a 0 to 100 scale
`ncuts = 100`. A scale that starts at 1 must be shifted to start at 0: a
Likert item 1 to 5 becomes 0 to 4, with `ncuts = 4`.

``` r

likert <- c(1, 2, 3, 5, 4, 1)
brs_check(likert - 1, ncuts = 4)
#>         left   right      yt y delta
#> [1,] 0.00001 0.12500 0.00001 0     1
#> [2,] 0.12500 0.37500 0.25000 1     3
#> [3,] 0.37500 0.62500 0.50000 2     3
#> [4,] 0.87500 0.99999 0.99999 4     2
#> [5,] 0.62500 0.87500 0.75000 3     3
#> [6,] 0.00001 0.12500 0.00001 0     1
```

Without the shift, category 1 would be read as an interior score and the
lower border of the scale would be a category that nobody can choose.

## The latent beta model

Each score is the coarsened reading of a latent $`Y \in (0, 1)`$ with a
beta law. The default parameterization (`repar = 2`) uses the mean
$`\mu`$ and the dispersion $`\sigma \in (0, 1)`$:
``` math
a = \mu\,\frac{1 - \sigma}{\sigma}, \qquad
b = (1 - \mu)\,\frac{1 - \sigma}{\sigma}, \qquad
\sigma = \frac{1}{1 + a + b} = \frac{1}{1 + \phi},
```
where $`\phi = a + b`$ is the usual beta precision. Then
``` math
E(Y) = \mu, \qquad \mathrm{Var}(Y) = \mu(1 - \mu)\,\sigma,
```
so $`\sigma`$ is the variance relative to its maximum $`\mu(1 - \mu)`$.
Both parameters can depend on covariates:
``` math
g(\mu_i) = x_i^\top \beta, \qquad h(\sigma_i) = z_i^\top \gamma .
```
The link `link` for $`\mu`$ and `link_phi` for $`\sigma`$ can be
`"logit"` (default), `"probit"`, `"cauchit"` or `"cloglog"`. A one-part
formula `y ~ x1 + x2` keeps $`\sigma`$ constant; `y ~ x1 + x2 | z1`
models it (see
[`vignette("brs-advanced-workflows")`](https://evandeilton.github.io/betaregscale/articles/brs-advanced-workflows.md)).

[`brs_repar()`](https://evandeilton.github.io/betaregscale/reference/brs_repar.md)
gives the shapes. In the code the second parameter is always called
`phi`; under `repar = 2` it is the dispersion $`\sigma`$:

``` r

brs_repar(mu = 0.3, phi = 0.2, repar = 2)         # mean 0.3, dispersion 0.2
#>   shape1 shape2
#> 1    1.2    2.8
brs_repar(mu = 0.3, phi = 1 / 0.2 - 1, repar = 1) # same law: precision 4
#>   shape1 shape2
#> 1    1.2    2.8
```

Two other parameterizations are available. Their admissible links
differ, and the first link of each cell is the default:

| `repar` | parameters | shapes $`(a, b)`$ | `link` | `link_phi` |
|:---|:---|:---|:---|:---|
| 0 | shapes $`p, q > 0`$ | $`(p, q)`$ | log, sqrt | log, sqrt |
| 1 | mean $`\mu`$, precision $`\phi > 0`$ | $`(\mu\phi, (1-\mu)\phi)`$ | logit, probit, cauchit, cloglog | log, sqrt |
| 2 | mean $`\mu`$, dispersion $`\sigma \in (0, 1)`$ | $`(\mu(1-\sigma)/\sigma, (1-\mu)(1-\sigma)/\sigma)`$ | logit, probit, cauchit, cloglog | logit, probit, cauchit, cloglog |

Any other combination is an error. For numerical safety the mean and the
dispersion are clamped to $`[10^{-5}, 1 - 10^{-5}]`$ inside the
likelihood.

## Interval maps and censoring types

A score $`s`$ is mapped to a cell of $`(0, 1)`$. The map is set by
`interval`:

| `interval` | cell of score $`s`$ | latent score of $`y^* \in (0, 1)`$ |
|:---|:---|:---|
| `"mid"` (default) | $`[s - h,\, s + h]/K`$ | $`K y^*`$ |
| `"right"` | $`[s,\, s + 1]/(K + 1)`$ | $`(K + 1) y^*`$ |
| `"left"` | $`[s,\, s + 1]/(K + 1)`$ | $`(K + 1) y^* - 1`$ |

$`h`$ is the half-width `lim` (default 0.5, used by `"mid"` only; values
below 0.5 leave gaps between cells and give a warning). Under `"mid"`
the cell of a score is centered on its label; under `"right"` and
`"left"` the $`K + 1`$ cells have equal width. `"right"` and `"left"`
give the same likelihood and differ only in the latent score returned by
`predict(type = "score")`.

[`brs_check()`](https://evandeilton.github.io/betaregscale/reference/brs_check.md)
returns the cell and the censoring type $`\delta`$ of each score: score
0 is left-censored ($`\delta = 1`$), score $`K`$ right-censored
($`\delta = 2`$), the others interval-censored ($`\delta = 3`$), and a
value already in $`(0, 1)`$ is an exact observation ($`\delta = 0`$):

``` r

brs_check(c(0, 3, 5, 10), ncuts = 10)
#>         left   right      yt  y delta
#> [1,] 0.00001 0.05000 0.00001  0     1
#> [2,] 0.25000 0.35000 0.30000  3     3
#> [3,] 0.45000 0.55000 0.50000  5     3
#> [4,] 0.95000 0.99999 0.99999 10     2
brs_check(c(0, 3, 5, 10), ncuts = 10, interval = "right")
#>           left      right         yt  y delta
#> [1,] 0.0000100 0.09090909 0.04545455  0     1
#> [2,] 0.2727273 0.36363636 0.31818182  3     3
#> [3,] 0.4545455 0.54545455 0.50000000  5     3
#> [4,] 0.9090909 0.99999000 0.95454545 10     2
```

The endpoints are clamped to $`[\epsilon, 1 - \epsilon]`$ with
$`\epsilon = 10^{-5}`$ only as a numerical guard, so the border scores
need no transformation that depends on the sample size. The cells are
fixed by the scale and not by the observed value, so the censoring is
non-informative.

## The likelihood is a pmf over the scores

With beta density $`f`$ and distribution function $`F`$ at the shapes
$`(a_i, b_i)`$ and cell $`[l_i, u_i]`$, observation $`i`$ contributes

| $`\delta`$ | observation                    | contribution $`L_i`$ |
|:----------:|:-------------------------------|:---------------------|
|     0      | exact value $`y_i \in (0, 1)`$ | $`f(y_i)`$           |
|     1      | score 0                        | $`F(u_i)`$           |
|     2      | score $`K`$                    | $`1 - F(l_i)`$       |
|     3      | interior score                 | $`F(u_i) - F(l_i)`$  |

and $`L(\beta, \gamma) = \prod_i L_i`$. For score data with `lim = 0.5`
the cells partition $`(0, 1)`$, so the contributions are the
probabilities $`P(S_i = s_i)`$ of a probability mass function over
$`0, \ldots, K`$. The log-likelihood is the log probability of the
observed scores, on the same scale as that of any other model for the
same scores.

A first fit on simulated NRS-11 scores (`phi` in
[`brs_sim()`](https://evandeilton.github.io/betaregscale/reference/brs_sim.md)
is given on the link scale, here $`\sigma = 0.15`$):

``` r

set.seed(4255)
d <- data.frame(x1 = rnorm(200), x2 = rbinom(200, 1, 0.5))
sim <- brs_sim(~ x1 + x2, data = d, beta = c(-0.4, 0.5, -0.6),
               phi = qlogis(0.15), ncuts = 10)
head(sim)
#>      left right    yt y delta          x1 x2
#> 1 0.65000  0.75 7e-01 7     3  1.95102377  0
#> 2 0.65000  0.75 7e-01 7     3  0.77253858  0
#> 3 0.35000  0.45 4e-01 4     3  0.72640816  0
#> 4 0.15000  0.25 2e-01 2     3  0.04873961  1
#> 5 0.00001  0.05 1e-05 0     1 -0.54450108  1
#> 6 0.45000  0.55 5e-01 5     3  0.36002855  0

fit <- brs(y ~ x1 + x2, data = sim)
summary(fit)
#> 
#> Call:
#> brs(formula = y ~ x1 + x2, data = sim)
#> 
#> Quantile residuals:
#>     Min      1Q  Median      3Q     Max 
#> -2.6209 -0.6469  0.0334  0.6900  2.9959 
#> 
#> Coefficients (mean model with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -0.34786    0.07834  -4.441 8.97e-06 ***
#> x1           0.50423    0.06043   8.343  < 2e-16 ***
#> x2          -0.69910    0.11277  -6.199 5.67e-10 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Coefficients (dispersion model for sigma with logit link):
#>       Estimate Std. Error z value Pr(>|z|)    
#> (phi)  -1.8827     0.1024  -18.38   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Log-likelihood: -371.6645 on 4 Df | AIC: 751.3290 | BIC: 764.5223 
#> Pseudo R-squared: 0.2721  (midpoint approx.; interpret with caution for heavily censored data) 
#> Number of iterations: 29 (BFGS) 
#> Censoring: 188 interval | 11 left | 1 right
```

The score probabilities of each row sum to 1, and the per-observation
log probabilities add up to
[`logLik()`](https://rdrr.io/r/stats/logLik.html):

``` r

P <- brs_predict_scoreprob(fit)
round(P[1:3, ], 3)
#>      score_0 score_1 score_2 score_3 score_4 score_5 score_6 score_7 score_8
#> [1,]   0.000   0.002   0.013   0.039   0.083   0.136   0.186   0.213   0.196
#> [2,]   0.001   0.018   0.065   0.124   0.173   0.195   0.183   0.139   0.078
#> [3,]   0.001   0.019   0.069   0.128   0.176   0.195   0.180   0.135   0.074
#>      score_9 score_10
#> [1,]   0.119    0.014
#> [2,]   0.024    0.001
#> [3,]   0.022    0.001
range(rowSums(P))
#> [1] 1 1
c(sum = sum(brs_loglik_contrib(fit)), logLik = as.numeric(logLik(fit)))
#>       sum    logLik 
#> -371.6645 -371.6645
```

## Estimation

[`optim()`](https://rdrr.io/r/stats/optim.html) (BFGS by default,
`maxit = 5000` unless `control` says otherwise) maximizes the
log-likelihood, evaluated in C++. The gradient and the default Hessian
(`hessian_method = "cpp"`) use the chain rule on the linear predictors
with per-observation numerical derivatives; standard errors come from
the inverse of the negative Hessian.
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html) and
[`confint()`](https://rdrr.io/r/stats/confint.html) take
`model = "mean"` or `"precision"` to select a submodel:

``` r

coef(fit, model = "mean")
#> (Intercept)          x1          x2 
#>  -0.3478562   0.5042255  -0.6990964
confint(fit)
#>                  2.5 %     97.5 %
#> (Intercept) -0.5013907 -0.1943216
#> x1           0.3857768  0.6226742
#> x2          -0.9201170 -0.4780758
#> (phi)       -2.0835272 -1.6819461
```

After the fit,
[`brs()`](https://evandeilton.github.io/betaregscale/reference/brs.md)
checks the gradient, the Hessian and the clamps, stores the numbers in
`fit$diagnostics` and gives a one-line warning for each problem:

| Message | Meaning | What to do |
|:---|:---|:---|
| model matrix is rank deficient (error) | aliased columns | remove the named columns |
| model matrix is nearly collinear | condition number above $`10^4`$ | drop, combine or center the named columns |
| Gradient not ~0 at the estimate | a Newton step would still gain more than 0.01 in log-likelihood | refit with the other `method` or rescaled covariates |
| Hessian near-singular or not negative definite | a direction is nearly flat | [`vcov()`](https://rdrr.io/r/stats/vcov.html) gives `NA` there; use likelihood-ratio tests or the bootstrap |
| observations on the clamp boundary | means or shapes at the clamps, e.g. a group entirely at one border | merge sparse groups or drop the separating covariate |

``` r

str(fit$diagnostics[c("grad_gain", "min_eig_scaled", "n_clamped")])
#> List of 3
#>  $ grad_gain     : num 1.19e-09
#>  $ min_eig_scaled: num 0.305
#>  $ n_clamped     : int 0
```

## Preparing data with `brs_prep()`

[`brs()`](https://evandeilton.github.io/betaregscale/reference/brs.md)
applies
[`brs_check()`](https://evandeilton.github.io/betaregscale/reference/brs_check.md)
to raw scores.
[`brs_prep()`](https://evandeilton.github.io/betaregscale/reference/brs_prep.md)
accepts, row by row, a score only, a score and a censoring type `delta`,
interval bounds `left`/`right` (one of them `NA` for one-sided
censoring) or a score with both bounds. Bounds are on the latent score
scale:

``` r

d_prep <- data.frame(
  y     = c(5, 0, NA, NA, 5),
  delta = c(3, NA, NA, NA, 0),
  left  = c(NA, NA, 2, NA, NA),
  right = c(NA, NA, 5, 3, NA)
)
brs_prep(d_prep, ncuts = 10)
#> brs_prep: n = 5 | exact = 1, left = 2, right = 0, interval = 2
#>      left right      yt   y delta
#> 1 0.45000  0.55 0.50000 5.0     3
#> 2 0.00001  0.05 0.00001 0.0     1
#> 3 0.20000  0.50 0.35000 3.5     3
#> 4 0.00001  0.30 0.15000 1.5     1
#> 5 0.50000  0.50 0.50000 5.0     0
```

The rows are, in order: an interior score, the lowest score, an analyst
interval $`[2, 5]`$, a bound “below 3” (left-censored) and a score read
as an exact value. When the data have `left`/`right` columns, an
interior score whose `left`, `right` and `delta` are all missing is read
as an exact value at the cell center; without those columns it is
interval-censored.
[`brs()`](https://evandeilton.github.io/betaregscale/reference/brs.md)
and
[`brsmm()`](https://evandeilton.github.io/betaregscale/reference/brsmm.md)
use prepared data as they are, with the stored `ncuts`, `lim` and
`interval`; so does the output of
[`brs_sim()`](https://evandeilton.github.io/betaregscale/reference/brs_sim.md).

## Choosing the interval map

The three maps are different coarsening models of the same scores. With
`lim = 0.5` each one gives a pmf over $`0, \ldots, K`$, so
log-likelihoods, AIC and BIC of fits under different maps are
comparable. The maps are not nested, so there is no likelihood-ratio
test between them and [`anova()`](https://rdrr.io/r/stats/anova.html)
stops:

``` r

raw <- data.frame(y = sim$y, x1 = sim$x1, x2 = sim$x2)
fits <- lapply(c(mid = "mid", right = "right", left = "left"), function(m)
  brs(y ~ x1 + x2, data = raw, ncuts = 10, interval = m))
data.frame(logLik = sapply(fits, logLik), AIC = sapply(fits, AIC),
           BIC = sapply(fits, BIC))
#>          logLik      AIC      BIC
#> mid   -371.6645 751.3290 764.5223
#> right -373.3736 754.7472 767.9405
#> left  -373.3736 754.7472 767.9405
tryCatch(anova(fits$mid, fits$right), error = conditionMessage)
#> [1] "All models must use the same `interval` (found: mid, right); fits with different interval maps are not nested; compare them with AIC() or BIC()."
```

The scores were simulated under `"mid"`, which has the largest
log-likelihood here; `"right"` and `"left"` give the same value.

## Boundary cutpoints

Under a fixed map, score $`s`$ has the cell $`(c_{s-1}, c_s)`$ with
cutpoints set by the score labels. Rating data often pile up at the ends
of the scale more than these cells allow. `cutpoints = "boundary"`
estimates the two border cutpoints, $`c_{lo}`$ between scores 0 and 1
and $`c_{hi}`$ between $`K - 1`$ and $`K`$, and keeps all the others.
With $`\epsilon = 10^{-5}`$ (the clamp above), $`u_1`$ the upper end of
the cell of score 1 and $`l_{K-1}`$ the lower end of the cell of score
$`K - 1`$,
``` math
c_{lo} = \epsilon + (u_1 - \epsilon)\,\mathrm{logit}^{-1}(g_{lo}), \qquad
c_{hi} = l_{K-1} + (1 - \epsilon - l_{K-1})\,\mathrm{logit}^{-1}(g_{hi}),
```
so the free parameters $`(g_{lo}, g_{hi})`$ are unconstrained and the
cutpoints stay ordered. The fixed map is an interior point of this
model, so `anova(fixed, boundary)` is a likelihood-ratio test on 2
degrees of freedom. It needs score data (no exact values or analyst
bounds), `lim = 0.5`, the scores 0, 1, $`K - 1`$ and $`K`$ observed, and
$`K \ge 3`$: with $`K = 2`$ the cutpoints are not identified and
[`brs()`](https://evandeilton.github.io/betaregscale/reference/brs.md)
stops, and with $`K = 3`$ they are weakly identified and a warning says
so. Use the test for $`K \ge 4`$.

Scores with extra mass at both ends, simulated with
`brs_sim(cutpoints =)`:

``` r

set.seed(1)
heap <- brs_sim(~ x1 + x2, data = d, beta = c(0, 0.5, -0.6),
                phi = qlogis(0.2), ncuts = 10, cutpoints = c(0.12, 0.88))
table(heap$y)
#> 
#>  0  1  2  3  4  5  6  7  8  9 10 
#> 29  8 18 33 30 24 19 17 14  3  5

fit_fix <- brs(y ~ x1 + x2, data = heap)
fit_bnd <- brs(y ~ x1 + x2, data = heap, cutpoints = "boundary")
anova(fit_fix, fit_bnd)
#> Likelihood-ratio comparison of brs/brsmm models
#> Rows M2: cutpoints "fixed" vs "boundary" (2 cutpoint parameters; the fixed map is an interior point; chi2 reference when both cutpoints are identified).
#> 
#>          Df  logLik    AIC    BIC  Chisq Chi Df Pr(>Chisq)    
#> M1 (brs)  4 -432.49 872.98 886.18                             
#> M2 (brs)  6 -415.62 843.24 863.03 33.738      2  4.719e-08 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
coef(fit_bnd, model = "cutpoints")
#>  (cut_lo)  (cut_hi) 
#>  1.401056 -1.023602
```

`summary(fit_bnd)` reports the cutpoints on $`(0, 1)`$ and on the score
scale with delta-method standard errors. Score probabilities, expected
scores, quantile residuals, the bootstrap and cross-validation use the
estimated cutpoints. A longer example is in
[`vignette("brs-advanced-workflows")`](https://evandeilton.github.io/betaregscale/articles/brs-advanced-workflows.md).

## References

- Lopes, J. E. (2023). *Modelos de regressão beta para dados de escala*.
  Master’s dissertation, Universidade Federal do Paraná, Curitiba.
  <https://hdl.handle.net/1884/86624>
- Bayer, F. M. (2011). *Modelagem e inferência em regressão beta*. PhD
  thesis, Universidade Federal de Pernambuco.
- Ferrari, S. L. P., and Cribari-Neto, F. (2004). Beta regression for
  modelling rates and proportions. *Journal of Applied Statistics*,
  31(7), 799–815. <https://doi.org/10.1080/0266476042000214501>
- Smithson, M., and Verkuilen, J. (2006). A better lemon squeezer?
  Maximum-likelihood regression with beta-distributed dependent
  variables. *Psychological Methods*, 11(1), 54–71.
  <https://doi.org/10.1037/1082-989X.11.1.54>
