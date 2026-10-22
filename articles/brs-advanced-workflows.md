# Advanced workflows: dispersion, cutpoints and model comparison

This vignette works through three modeling questions: whether the
dispersion depends on covariates, whether the border scores need their
own cutpoints, and how to compare models that differ in link, interval
map or model family. The model and the notation are those of
[`vignette("brs-intro")`](https://evandeilton.github.io/betaregscale/articles/brs-intro.md).
The data are simulated and small, so the results are illustrative.

``` r

library(betaregscale)
```

## Variable dispersion

A two-part formula `y ~ x1 + x2 | z1` models the dispersion with its own
linear predictor, $`h(\sigma_i) = z_i^\top \gamma`$ (logit link by
default under `repar = 2`).
[`brs_sim()`](https://evandeilton.github.io/betaregscale/reference/brs_sim.md)
takes the dispersion coefficients in `zeta`, on the link scale:

``` r

set.seed(11)
d <- data.frame(x1 = rnorm(300), x2 = rbinom(300, 1, 0.5), z1 = rnorm(300))
sv <- brs_sim(~ x1 + x2 | z1, data = d, beta = c(-0.2, 0.5, -0.5),
              zeta = c(qlogis(0.15), 0.6), ncuts = 10)
table(sv$y)
#> 
#>  0  1  2  3  4  5  6  7  8  9 10 
#> 11 35 40 55 31 42 35 24 13 12  2
```

The constant-dispersion model is nested in the variable one, so
[`anova()`](https://rdrr.io/r/stats/anova.html) gives a likelihood-ratio
test on the added dispersion coefficients:

``` r

fit_const <- brs(y ~ x1 + x2, data = sv)
fit_var <- brs(y ~ x1 + x2 | z1, data = sv)
anova(fit_const, fit_var)
#> Likelihood-ratio comparison of brs/brsmm models
#> 
#>          Df  logLik    AIC    BIC  Chisq Chi Df Pr(>Chisq)    
#> M1 (brs)  4 -611.69 1231.4 1246.2                             
#> M2 (brs)  5 -588.92 1187.8 1206.4 45.541      1  1.495e-11 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
summary(fit_var)
#> 
#> Call:
#> brs(formula = y ~ x1 + x2 | z1, data = sv)
#> 
#> Quantile residuals:
#>     Min      1Q  Median      3Q     Max 
#> -2.7415 -0.6239  0.0050  0.6938  2.8562 
#> 
#> Coefficients (mean model with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -0.16331    0.06319  -2.585  0.00975 ** 
#> x1           0.55034    0.04889  11.256  < 2e-16 ***
#> x2          -0.47818    0.08873  -5.389 7.08e-08 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Coefficients (dispersion model for sigma with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -1.79655    0.08254 -21.767  < 2e-16 ***
#> z1           0.60411    0.09071   6.659 2.75e-11 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Log-likelihood: -588.9232 on 5 Df | AIC: 1187.8464 | BIC: 1206.3653 
#> Pseudo R-squared: 0.1761  (midpoint approx.; interpret with caution for heavily censored data) 
#> Number of iterations: 39 (BFGS) 
#> Censoring: 287 interval | 11 left | 2 right
```

The dispersion coefficients carry the `(phi)_` prefix;
`coef(fit, model = "precision")` selects them, and
`predict(type = "precision")` returns $`\sigma_i`$:

``` r

coef(fit_var, model = "precision")
#> (phi)_(Intercept)          (phi)_z1 
#>        -1.7965462         0.6041058
nd <- data.frame(x1 = 0, x2 = 0, z1 = c(-1, 0, 1))
cbind(nd, sigma = predict(fit_var, nd, type = "precision"),
      sd_latent = 10 * sqrt(predict(fit_var, nd, type = "variance")))
#>   x1 x2 z1      sigma sd_latent
#> 1  0  0 -1 0.08312299  1.436760
#> 2  0  0  0 0.14227201  1.879678
#> 3  0  0  1 0.23282275  2.404564
```

`brs_marginaleffects(fit, model = "precision")` gives average marginal
effects on the dispersion submodel.

## Boundary cutpoints

Rating data often pile up at 0 and at $`K`$. With
`cutpoints = "boundary"` the cutpoints between scores 0 and 1 and
between $`K - 1`$ and $`K`$ are estimated (see
[`vignette("brs-intro")`](https://evandeilton.github.io/betaregscale/articles/brs-intro.md)
for the parameterization and the requirements). Scores simulated with
border cutpoints 0.12 and 0.88 instead of the 0.05 and 0.95 of the
`"mid"` map:

``` r

set.seed(1)
heap <- brs_sim(~ x1 + x2, data = d, beta = c(0, 0.5, -0.6),
                phi = qlogis(0.2), ncuts = 10, cutpoints = c(0.12, 0.88))
table(heap$y)
#> 
#>  0  1  2  3  4  5  6  7  8  9 10 
#> 33  9 32 38 44 36 40 27 24  5 12
```

The fixed map is nested at an interior point of the boundary model, so
the likelihood-ratio test has 2 degrees of freedom ($`K \ge 4`$):

``` r

fit_fix <- brs(y ~ x1 + x2, data = heap)
fit_bnd <- brs(y ~ x1 + x2, data = heap, cutpoints = "boundary")
anova(fit_fix, fit_bnd)
#> Likelihood-ratio comparison of brs/brsmm models
#> Rows M2: cutpoints "fixed" vs "boundary" (2 cutpoint parameters; the fixed map is an interior point; chi2 reference when both cutpoints are identified).
#> 
#>          Df  logLik    AIC    BIC  Chisq Chi Df Pr(>Chisq)    
#> M1 (brs)  4 -664.39 1336.8 1351.6                             
#> M2 (brs)  6 -636.12 1284.2 1306.5 56.541      2  5.276e-13 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

[`summary()`](https://rdrr.io/r/base/summary.html) reports the estimated
cutpoints on $`(0, 1)`$ and on the score scale, with delta-method
standard errors, next to their fixed values;
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html) and
[`confint()`](https://rdrr.io/r/stats/confint.html) accept
`model = "cutpoints"`:

``` r

summary(fit_bnd)
#> 
#> Call:
#> brs(formula = y ~ x1 + x2, data = heap, cutpoints = "boundary")
#> 
#> Quantile residuals:
#>     Min      1Q  Median      3Q     Max 
#> -2.4684 -0.6102 -0.0448  0.6593  2.7993 
#> 
#> Coefficients (mean model with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.02673    0.07589   0.352    0.725    
#> x1           0.49265    0.05868   8.395  < 2e-16 ***
#> x2          -0.47655    0.10629  -4.483 7.35e-06 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Coefficients (dispersion model for sigma with logit link):
#>       Estimate Std. Error z value Pr(>|z|)    
#> (phi) -1.42562    0.08684  -16.42   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Boundary cutpoints (logit scale; on (0, 1) and on the score scale, with delta-method SEs; 'fixed' = the map's value):
#>             Estimate Std. Error cutpoint SE(cut)  score SE(score) fixed
#> c_lo (0|1)    1.5165     0.3784   0.1230  0.0084 1.2301    0.0838   0.5
#> c_hi (9|10)  -1.3994     0.5067   0.8797  0.0121 8.7969    0.1206   9.5
#> 
#> ---
#> Log-likelihood: -636.1209 on 6 Df | AIC: 1284.2418 | BIC: 1306.4645 
#> Pseudo R-squared: 0.2617  (midpoint approx.; interpret with caution for heavily censored data) 
#> Number of iterations: 35 (BFGS) 
#> Censoring: 255 interval | 33 left | 12 right
confint(fit_bnd, model = "cutpoints")
#>               2.5 %     97.5 %
#> (cut_lo)  0.7748769  2.2582010
#> (cut_hi) -2.3924205 -0.4063316
```

The average predicted score distribution against the observed one shows
what the two cutpoints change:

``` r

obs <- prop.table(table(factor(heap$y, levels = 0:10)))
round(cbind(observed = as.numeric(obs),
            fixed = colMeans(brs_predict_scoreprob(fit_fix)),
            boundary = colMeans(brs_predict_scoreprob(fit_bnd))), 3)
#>          observed fixed boundary
#> score_0     0.110 0.070    0.105
#> score_1     0.030 0.121    0.030
#> score_2     0.107 0.118    0.123
#> score_3     0.127 0.115    0.132
#> score_4     0.147 0.111    0.133
#> score_5     0.120 0.106    0.127
#> score_6     0.133 0.099    0.115
#> score_7     0.090 0.090    0.098
#> score_8     0.080 0.079    0.076
#> score_9     0.017 0.065    0.018
#> score_10    0.040 0.026    0.043
```

Cross-validation refits the cutpoints in each training fold and scores
the held-out rows with them. A fold without one of the scores 0, 1,
$`K - 1`$ or $`K`$ cannot be fitted and is reported in the `error`
column.

``` r

set.seed(2)
cv_fix <- brs_cv(y ~ x1 + x2, data = heap, k = 5)
set.seed(2)
cv_bnd <- brs_cv(y ~ x1 + x2, data = heap, k = 5, cutpoints = "boundary")
c(fixed = mean(cv_fix$log_score), boundary = mean(cv_bnd$log_score))
#>     fixed  boundary 
#> -2.222953 -2.138410
```

Pseudo $`R^2`$ and the midpoint residuals use the centers of the
estimated border cells, so they are not directly comparable between
fixed-map and boundary fits; compare those fits by likelihood or
held-out log score.

## Model comparison

### Links and interval maps

Fits with different links for the same scores are not nested; compare
them by AIC or BIC:

``` r

fits <- lapply(c(logit = "logit", probit = "probit", cloglog = "cloglog"),
               function(l) brs(y ~ x1 + x2 | z1, data = sv, link = l))
brs_table(models = fits, sort_by = "AIC")
#>     model nobs npar    logLik      AIC      BIC pseudo_r2 exact left right
#> 1 cloglog  300    5 -587.8171 1185.634 1204.153    0.1454     0   11     2
#> 2   logit  300    5 -588.9232 1187.846 1206.365    0.1761     0   11     2
#> 3  probit  300    5 -589.1080 1188.216 1206.735    0.2399     0   11     2
#>   interval prop_exact prop_left prop_right prop_interval
#> 1      287          0    0.0367     0.0067        0.9567
#> 2      287          0    0.0367     0.0067        0.9567
#> 3      287          0    0.0367     0.0067        0.9567
```

The same holds across interval maps: with `lim = 0.5` each map gives a
pmf over $`0, \ldots, K`$, so
[`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html) are comparable across
`interval = "mid"`, `"right"` and `"left"`, but the maps are not nested
and [`anova()`](https://rdrr.io/r/stats/anova.html) refuses a
likelihood-ratio test between them (see
[`vignette("brs-intro")`](https://evandeilton.github.io/betaregscale/articles/brs-intro.md)).

### Comparison with an ordinal model

Because the likelihood of
[`brs()`](https://evandeilton.github.io/betaregscale/reference/brs.md)
is the probability of the observed scores, it is on the same scale as
the likelihood of an ordinal regression for the same scores, and the two
can be compared by AIC. A cumulative logit model
$`P(S \le s) = \mathrm{logit}^{-1}(\theta_s - x^\top
\beta)`$, with $`K`$ free thresholds, can be fitted in a few lines of
base R:

``` r

X <- model.matrix(~ x1 + x2, data = heap)[, -1]
s <- heap$y
K <- 10
# Negative log-likelihood; thresholds kept ordered by positive increments
nll_ord <- function(par) {
  theta <- cumsum(c(par[1], exp(par[2:K])))
  eta <- drop(X %*% par[-(1:K)])
  cum <- cbind(0, plogis(outer(-eta, theta, "+")), 1)
  -sum(log(cum[cbind(seq_along(s), s + 2)] - cum[cbind(seq_along(s), s + 1)]))
}
start <- c(qlogis(mean(s == 0)), rep(log(0.4), K - 1), 0, 0)
ord <- optim(start, nll_ord, method = "BFGS", control = list(maxit = 1000))
ord$convergence
#> [1] 0

ll <- c(brs_fixed = as.numeric(logLik(fit_fix)),
        brs_boundary = as.numeric(logLik(fit_bnd)),
        cumulative_logit = -ord$value)
df <- c(length(coef(fit_fix)), length(coef(fit_bnd)), length(ord$par))
data.frame(logLik = ll, df = df, AIC = -2 * ll + 2 * df)
#>                     logLik df      AIC
#> brs_fixed        -664.3913  4 1336.783
#> brs_boundary     -636.1209  6 1284.242
#> cumulative_logit -636.4632 12 1296.926
```

The cumulative logit model spends $`K`$ parameters on the thresholds;
the beta model with boundary cutpoints uses two plus the dispersion. The
same reasoning applies out of sample: `brs_loglik_contrib(fit, newdata)`
and
[`brsmm_cv()`](https://evandeilton.github.io/betaregscale/reference/brsmm_cv.md)
return held-out log probabilities of the scores that can be set against
those of any model giving a pmf on $`0, \ldots, K`$, evaluated on the
same rows or folds.

## A Monte Carlo study

A simulation study of the estimator is a loop of
[`brs_sim()`](https://evandeilton.github.io/betaregscale/reference/brs_sim.md)
and
[`brs()`](https://evandeilton.github.io/betaregscale/reference/brs.md)
over a fixed design. The code below is a template (not run here): logit
links, `interval = "mid"`, two scale lengths, two sample sizes and three
dispersion levels, with bias, root mean squared error and coverage of
the 95% Wald interval per cell. Increase `R` for precise coverages; with
`R = 200` the Monte Carlo standard error of a coverage near 0.95 is
about 0.015.

``` r

R <- 200L
beta <- c(0.5, 0.2, -0.7)
grid <- expand.grid(sigma = c(0.08, 0.2, 0.4), n = c(100L, 500L),
                    K = c(10L, 100L))
one_cell <- function(sigma, n, K) {
  set.seed(1000 + n + K + round(100 * sigma))
  d <- data.frame(x1 = rnorm(n), x2 = rbinom(n, 1, 0.5))   # fixed design
  truth <- c(beta, qlogis(sigma))
  est <- se <- matrix(NA_real_, R, 4)
  for (r in seq_len(R)) {
    s <- brs_sim(~ x1 + x2, data = d, beta = beta, phi = qlogis(sigma),
                 ncuts = K)
    f <- tryCatch(brs(y ~ x1 + x2, data = s), error = function(e) NULL)
    if (is.null(f) || f$convergence != 0L) next
    est[r, ] <- coef(f)
    se[r, ] <- sqrt(diag(vcov(f)))
  }
  ok <- stats::complete.cases(est, se)
  err <- est[ok, ] - rep(truth, each = sum(ok))
  data.frame(sigma, n, K, par = c("b0", "b1", "b2", "logit_sigma"),
             bias = colMeans(err), rmse = sqrt(colMeans(err^2)),
             cover = colMeans(abs(err) <= qnorm(0.975) * se[ok, ]))
}
res <- do.call(rbind, Map(one_cell, grid$sigma, grid$n, grid$K))
```

For mixed models the same loop uses
[`brsmm_sim()`](https://evandeilton.github.io/betaregscale/reference/brsmm_sim.md)
and
[`brsmm()`](https://evandeilton.github.io/betaregscale/reference/brsmm.md);
the attribute `attr(sim, "truth")$par` holds the true parameters in the
order of [`coef()`](https://rdrr.io/r/stats/coef.html).

## References

- Lopes, J. E. (2023). *Modelos de regressão beta para dados de escala*.
  Master’s dissertation, Universidade Federal do Paraná, Curitiba.
  <https://hdl.handle.net/1884/86624>
- McCullagh, P. (1980). Regression models for ordinal data. *Journal of
  the Royal Statistical Society: Series B*, 42(2), 109–142.
