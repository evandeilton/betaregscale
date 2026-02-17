# Advanced Workflows for High-Level Users

## Audience and scope

This vignette is designed for analysts who already know beta regression
and need a production-oriented workflow with:

1.  reproducible model selection;
2.  inferential and predictive diagnostics;
3.  mixed-effects escalation (`brs` -\> `brsmm`) with LR comparison;
4.  high-signal reporting tables for technical decision-making.

``` r
library(betaregscale)
```

## 1) Reproducible simulation and data checks

``` r
n <- 260
d <- data.frame(
  x1 = rnorm(n),
  x2 = rnorm(n),
  z1 = rnorm(n)
)

sim <- brs_sim(
  formula = ~ x1 + x2 | z1,
  data = d,
  beta = c(0.15, 0.55, -0.30),
  zeta = c(-0.10, 0.35),
  link = "logit",
  link_phi = "logit",
  ncuts = 100,
  repar = 2
)

knitr::kable(head(sim), digits = 4)
```

|  left | right |   yt |   y | delta |      x1 |      x2 |      z1 |
|------:|------:|-----:|----:|------:|--------:|--------:|--------:|
| 0.745 | 0.755 | 0.75 |  75 |     3 |  0.5206 |  0.2576 | -0.9272 |
| 0.705 | 0.715 | 0.71 |  71 |     3 | -1.0797 |  0.6131 |  0.7927 |
| 0.045 | 0.055 | 0.05 |   5 |     3 |  0.1392 | -0.6523 | -0.6891 |
| 0.025 | 0.035 | 0.03 |   3 |     3 | -0.0847 | -0.1383 | -0.0417 |
| 0.275 | 0.285 | 0.28 |  28 |     3 | -0.6666 | -0.1236 | -0.7856 |
| 0.000 | 0.005 | 0.00 |   0 |     1 | -2.5161 | -0.9647 |  0.4345 |

``` r
knitr::kable(
  data.frame(
    n = nrow(sim),
    exact = sum(sim$delta == 0),
    left = sum(sim$delta == 1),
    right = sum(sim$delta == 2),
    interval = sum(sim$delta == 3)
  ),
  digits = 4
)
```

|   n | exact | left | right | interval |
|----:|------:|-----:|------:|---------:|
| 260 |     0 |   15 |    16 |      229 |

## 2) Fixed-effects candidate set and model ranking

``` r
fit_logit <- brs(y ~ x1 + x2 | z1, data = sim, link = "logit", repar = 2)
fit_probit <- brs(y ~ x1 + x2 | z1, data = sim, link = "probit", repar = 2)
fit_cauchit <- brs(y ~ x1 + x2 | z1, data = sim, link = "cauchit", repar = 2)

tab_rank <- brs_table(
  logit = fit_logit,
  probit = fit_probit,
  cauchit = fit_cauchit
)
knitr::kable(tab_rank, digits = 4)
```

| model   | nobs | npar |    logLik |      AIC |      BIC | pseudo_r2 | exact | left | right | interval | prop_exact | prop_left | prop_right | prop_interval |
|:--------|-----:|-----:|----------:|---------:|---------:|----------:|------:|-----:|------:|---------:|-----------:|----------:|-----------:|--------------:|
| logit   |  260 |    5 | -1120.320 | 2250.639 | 2268.443 |    0.1875 |     0 |   15 |    16 |      229 |          0 |    0.0577 |     0.0615 |        0.8808 |
| probit  |  260 |    5 | -1120.281 | 2250.563 | 2268.366 |    0.1807 |     0 |   15 |    16 |      229 |          0 |    0.0577 |     0.0615 |        0.8808 |
| cauchit |  260 |    5 | -1120.656 | 2251.311 | 2269.115 |    0.1580 |     0 |   15 |    16 |      229 |          0 |    0.0577 |     0.0615 |        0.8808 |

## 3) Inference stack: Wald + bootstrap + AME

``` r
knitr::kable(brs_est(fit_logit), digits = 4)
```

| variable           | estimate |     se | z_value | p_value | ci_lower | ci_upper |
|:-------------------|---------:|-------:|--------:|--------:|---------:|---------:|
| (Intercept)        |   0.0673 | 0.0783 |  0.8589 |  0.3904 |  -0.0862 |   0.2208 |
| x1                 |   0.5759 | 0.0847 |  6.7979 |  0.0000 |   0.4099 |   0.7419 |
| x2                 |  -0.2344 | 0.0751 | -3.1193 |  0.0018 |  -0.3816 |  -0.0871 |
| (phi)\_(Intercept) |  -0.1677 | 0.0769 | -2.1800 |  0.0293 |  -0.3185 |  -0.0169 |
| (phi)\_z1          |   0.3492 | 0.0893 |  3.9085 |  0.0001 |   0.1741 |   0.5243 |

``` r
knitr::kable(confint(fit_logit), digits = 4)
```

|                    |   2.5 % |  97.5 % |
|:-------------------|--------:|--------:|
| (Intercept)        | -0.0862 |  0.2208 |
| x1                 |  0.4099 |  0.7419 |
| x2                 | -0.3816 | -0.0871 |
| (phi)\_(Intercept) | -0.3185 | -0.0169 |
| (phi)\_z1          |  0.1741 |  0.5243 |

``` r

boot_tab <- brs_bootstrap(fit_logit, R = 80, level = 0.95)
knitr::kable(head(boot_tab), digits = 4)
```

| parameter          | estimate | se_boot | ci_lower | ci_upper | level |
|:-------------------|---------:|--------:|---------:|---------:|------:|
| (Intercept)        |   0.0673 |  0.0765 |  -0.1117 |   0.1797 |  0.95 |
| x1                 |   0.5759 |  0.0786 |   0.4619 |   0.7419 |  0.95 |
| x2                 |  -0.2344 |  0.0736 |  -0.3600 |  -0.1112 |  0.95 |
| (phi)\_(Intercept) |  -0.1677 |  0.0699 |  -0.3160 |  -0.0475 |  0.95 |
| (phi)\_z1          |   0.3492 |  0.0805 |   0.2168 |   0.4992 |  0.95 |

``` r

ame_mu <- brs_marginaleffects(
  fit_logit,
  model = "mean",
  type = "response",
  interval = TRUE,
  n_sim = 120,
  seed = 2026
)
knitr::kable(ame_mu, digits = 4)
```

| variable |     ame | std.error | ci.lower | ci.upper | model | type     |   n |
|:---------|--------:|----------:|---------:|---------:|:------|:---------|----:|
| x1       |  0.1325 |    0.0170 |   0.1004 |   0.1652 | mean  | response | 260 |
| x2       | -0.0539 |    0.0154 |  -0.0804 |  -0.0232 | mean  | response | 260 |

## 4) Prediction layer on analyst scale

``` r
score_prob <- brs_predict_scoreprob(fit_logit, scores = 0:10)
knitr::kable(score_prob[1:8, 1:7], digits = 4)
```

| score_0 | score_1 | score_2 | score_3 | score_4 | score_5 | score_6 |
|--------:|--------:|--------:|--------:|--------:|--------:|--------:|
|  0.0048 |  0.0087 |  0.0084 |  0.0082 |  0.0081 |  0.0080 |  0.0079 |
|  0.1651 |  0.0641 |  0.0379 |  0.0284 |  0.0233 |  0.0200 |  0.0176 |
|  0.0068 |  0.0108 |  0.0099 |  0.0094 |  0.0090 |  0.0088 |  0.0086 |
|  0.0258 |  0.0250 |  0.0189 |  0.0162 |  0.0145 |  0.0134 |  0.0125 |
|  0.0267 |  0.0289 |  0.0226 |  0.0197 |  0.0179 |  0.0166 |  0.0157 |
|  0.2515 |  0.0772 |  0.0437 |  0.0319 |  0.0257 |  0.0217 |  0.0190 |
|  0.0354 |  0.0341 |  0.0256 |  0.0218 |  0.0195 |  0.0179 |  0.0167 |
|  0.1894 |  0.0709 |  0.0415 |  0.0310 |  0.0253 |  0.0216 |  0.0190 |

## 5) Out-of-sample validation

``` r
cv_tab <- brs_cv(
  y ~ x1 + x2 | z1,
  data = sim,
  k = 3,
  repeats = 1,
  seed = 2026,
  repar = 2
)
knitr::kable(cv_tab, digits = 4)
```

| repeat | fold | n_train | n_test | log_score | rmse_yt | mae_yt | converged | error |
|-------:|-----:|--------:|-------:|----------:|--------:|-------:|:----------|:------|
|      1 |    1 |     173 |     87 |   -4.5187 |  0.3555 | 0.2930 | TRUE      | NA    |
|      1 |    2 |     173 |     87 |   -4.3272 |  0.3326 | 0.2912 | TRUE      | NA    |
|      1 |    3 |     174 |     86 |   -4.1839 |  0.3073 | 0.2701 | TRUE      | NA    |

## 6) Escalation to mixed models

### 6.1 Simulate clustered data with random intercept + slope

``` r
g <- 16
ni <- 12
id <- factor(rep(seq_len(g), each = ni))
n_mm <- length(id)
x1 <- rnorm(n_mm)

b0 <- rnorm(g, sd = 0.40)
b1 <- rnorm(g, sd = 0.22)

eta_mu <- 0.20 + 0.65 * x1 + b0[id] + b1[id] * x1
eta_phi <- rep(-0.20, n_mm)

mu <- plogis(eta_mu)
phi <- plogis(eta_phi)
shp <- brs_repar(mu = mu, phi = phi, repar = 2)
y <- round(stats::rbeta(n_mm, shp$shape1, shp$shape2) * 100)

dmm <- data.frame(y = y, x1 = x1, id = id)
```

### 6.2 Fit evolutionary sequence

``` r
fit_brs <- brs(y ~ x1, data = dmm, repar = 2)
fit_ri <- brsmm(y ~ x1, random = ~ 1 | id, data = dmm, repar = 2)
fit_rs <- brsmm(y ~ x1, random = ~ 1 + x1 | id, data = dmm, repar = 2)

tab_lr <- anova(fit_brs, fit_ri, fit_rs, test = "Chisq")
knitr::kable(data.frame(model = rownames(tab_lr), tab_lr, row.names = NULL), digits = 4)
```

| model      |  Df |    logLik |      AIC |      BIC |   Chisq | Chi.Df | Pr..Chisq. |
|:-----------|----:|----------:|---------:|---------:|--------:|-------:|-----------:|
| M1 (brs)   |   3 | -829.5405 | 1665.081 | 1674.853 |      NA |     NA |         NA |
| M2 (brsmm) |   4 | -824.0537 | 1656.107 | 1669.137 | 10.9737 |      1 |     0.0009 |
| M3 (brsmm) |   6 | -821.7811 | 1655.562 | 1675.107 |  4.5452 |      2 |     0.1030 |

### 6.3 Model choice by LLR/LRT (ANOVA)

The [`anova()`](https://rdrr.io/r/stats/anova.html) methods provide a
practical LLR workflow:

- `M0 = brs` (no random effects);
- `M1 = brsmm` with random intercept (`~ 1 | id`);
- `M2 = brsmm` with random intercept + slope (`~ 1 + x1 | id`).

In nested comparisons, the statistic is:
$$LR = 2\{\ell\left( {\widehat{\theta}}_{\text{complex}} \right) - \ell\left( {\widehat{\theta}}_{\text{simple}} \right)\}.$$

For the first step (`M0 -> M1`), the null involves variance components
at the boundary ($\sigma_{b}^{2} = 0$); p-values should be interpreted
with caution. For `M1 -> M2`, the chi-square approximation is often used
as a practical decision aid.

``` r
tab_lr_df <- data.frame(model = rownames(tab_lr), tab_lr, row.names = NULL)
tab_lr_df$decision <- c(
  "baseline",
  ifelse(is.na(tab_lr_df$`Pr(>Chisq)`[2]), "inspect AIC/BIC + diagnostics",
         ifelse(tab_lr_df$`Pr(>Chisq)`[2] < 0.05, "prefer M1 over M0", "prefer M0 (parsimony)")),
  ifelse(is.na(tab_lr_df$`Pr(>Chisq)`[3]), "inspect AIC/BIC + diagnostics",
         ifelse(tab_lr_df$`Pr(>Chisq)`[3] < 0.05, "prefer M2 over M1", "prefer M1 (parsimony)"))
)
knitr::kable(tab_lr_df, digits = 4)
```

| model      |  Df |    logLik |      AIC |      BIC |   Chisq | Chi.Df | Pr..Chisq. | decision |
|:-----------|----:|----------:|---------:|---------:|--------:|-------:|-----------:|:---------|
| M1 (brs)   |   3 | -829.5405 | 1665.081 | 1674.853 |      NA |     NA |         NA | baseline |
| M2 (brsmm) |   4 | -824.0537 | 1656.107 | 1669.137 | 10.9737 |      1 |     0.0009 | baseline |
| M3 (brsmm) |   6 | -821.7811 | 1655.562 | 1675.107 |  4.5452 |      2 |     0.1030 | baseline |

## 7) Random-effects study (numeric + visual)

``` r
rs <- brsmm_re_study(fit_rs)
print(rs)
#> 
#> Random-effects study
#> Groups: 16 
#> 
#> Summary by term:
#>         term sd_model mean_mode sd_mode shrinkage_ratio shapiro_p
#>  (Intercept)   0.5127   -0.0038  0.4315          0.7083    0.8166
#>           x1   0.3745   -0.0114  0.2643          0.4982    0.9917
#> 
#> Estimated covariance matrix D:
#>        [,1]   [,2]
#> [1,] 0.2629 0.0369
#> [2,] 0.0369 0.1402
#> 
#> Estimated correlation matrix:
#>        [,1]   [,2]
#> [1,] 1.0000 0.1923
#> [2,] 0.1923 1.0000
knitr::kable(rs$summary, digits = 4)
```

| term        | sd_model | mean_mode | sd_mode | shrinkage_ratio | shapiro_p |
|:------------|---------:|----------:|--------:|----------------:|----------:|
| (Intercept) |   0.5127 |   -0.0038 |  0.4315 |          0.7083 |    0.8166 |
| x1          |   0.3745 |   -0.0114 |  0.2643 |          0.4982 |    0.9917 |

``` r
knitr::kable(rs$D, digits = 4)
```

|        |        |
|-------:|-------:|
| 0.2629 | 0.0369 |
| 0.0369 | 0.1402 |

``` r
knitr::kable(rs$Corr, digits = 4)
```

|        |        |
|-------:|-------:|
| 1.0000 | 0.1923 |
| 0.1923 | 1.0000 |

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  autoplot.brsmm(fit_rs, type = "ranef_caterpillar")
  autoplot.brsmm(fit_rs, type = "ranef_density")
  autoplot.brsmm(fit_rs, type = "ranef_pairs")
  autoplot.brsmm(fit_rs, type = "ranef_qq")
}
```

![](brs-advanced-workflows_files/figure-html/ranef-plots-1.png)

## 8) Practical decision checklist

- Start with `brs` candidates (`link`, `repar`) and rank with
  `brs_table`.
- Add bootstrap and AME before escalating complexity.
- Use `brs_cv` for out-of-sample stability checks.
- Escalate to `brsmm` only when LR/AIC/BIC and diagnostics support it.
- For random slopes, always inspect
  [`brsmm_re_study()`](https://evandeilton.github.io/betaregscale/reference/brsmm_re_study.md)
  and random-effects plots.

## References

Ferrari, S. L. P. and Cribari-Neto, F. (2004). Beta regression for
modelling rates and proportions. *Journal of Applied Statistics*, 31(7),
799-815. DOI: 10.1080/0266476042000214501. Validated online via:
<https://doi.org/10.1080/0266476042000214501>.

Smithson, M., and Verkuilen, J. (2006). A better lemon squeezer?
Maximum-likelihood regression with beta-distributed dependent variables.
*Psychological Methods*, 11(1), 54-71. DOI: 10.1037/1082-989X.11.1.54.
Validated online via: <https://doi.org/10.1037/1082-989X.11.1.54>.

Rue, H., Martino, S., and Chopin, N. (2009). Approximate Bayesian
inference for latent Gaussian models by using integrated nested Laplace
approximations. *Journal of the Royal Statistical Society: Series B*,
71(2), 319-392. DOI: 10.1111/j.1467-9868.2008.00700.x. Validated online
via: <https://doi.org/10.1111/j.1467-9868.2008.00700.x>.
