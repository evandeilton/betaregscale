# Mixed-Effects Beta Interval Regression with brsmm

## Overview

[`brsmm()`](https://evandeilton.github.io/betaregscale/reference/brsmm.md)
extends
[`brs()`](https://evandeilton.github.io/betaregscale/reference/brs.md)
to clustered data by adding Gaussian random effects in the mean submodel
while preserving the interval-censored beta likelihood for scale-derived
outcomes.

This vignette covers:

1.  full model mathematics;
2.  estimation by marginal maximum likelihood (Laplace approximation);
3.  practical use of all current `brsmm` methods;
4.  inferential and validation workflows, including parameter recovery.

``` r

library(betaregscale)
```

## Mathematical model

Assume observations $`i = 1, \dots, n_j`$ within groups
$`j = 1, \dots, G`$, with group-specific random-effects vector
$`\mathbf{b}_j \in \mathbb{R}^{q_b}`$.

### Linear predictors

``` math
\eta_{\mu,ij} = x_{ij}^\top \beta + w_{ij}^\top \mathbf{b}_j,
\qquad
\eta_{\phi,ij} = z_{ij}^\top \gamma
```

``` math
\mu_{ij} = g^{-1}(\eta_{\mu,ij}),
\qquad
\phi_{ij} = h^{-1}(\eta_{\phi,ij})
```

with $`g(\cdot)`$ and $`h(\cdot)`$ chosen by `link` and `link_phi`. The
random-effects design row $`w_{ij}`$ is defined by
`random = ~ terms | group`.

### Beta parameterization

For each $`(\mu_{ij},\phi_{ij})`$, `repar` maps to beta shape parameters
$`(a_{ij},b_{ij})`$ via
[`brs_repar()`](https://evandeilton.github.io/betaregscale/reference/brs_repar.md).

### Conditional contribution by censoring type

Each observation contributes:

``` math
L_{ij}(b_j;\theta)=
\begin{cases}
f(y_{ij}; a_{ij}, b_{ij}), & \delta_{ij}=0\\
F(u_{ij}; a_{ij}, b_{ij}), & \delta_{ij}=1\\
1 - F(l_{ij}; a_{ij}, b_{ij}), & \delta_{ij}=2\\
F(u_{ij}; a_{ij}, b_{ij}) - F(l_{ij}; a_{ij}, b_{ij}), & \delta_{ij}=3
\end{cases}
```

where $`l_{ij},u_{ij}`$ are interval endpoints on $`(0,1)`$,
$`f(\cdot)`$ is beta density, and $`F(\cdot)`$ is beta CDF.

### Random-effects distribution

``` math
\mathbf{b}_j \sim \mathcal{N}(\mathbf{0}, D),
```

where $`D`$ is a symmetric positive-definite covariance matrix.
Internally,
[`brsmm()`](https://evandeilton.github.io/betaregscale/reference/brsmm.md)
optimizes a packed lower-Cholesky parameterization $`D = LL^\top`$
(diagonal entries on log-scale for positivity).

### Group marginal likelihood

``` math
L_j(\theta)=\int_{\mathbb{R}^{q_b}}
\prod_{i=1}^{n_j} L_{ij}(b_j;\theta)\,
\varphi_{q_b}(\mathbf{b}_j;\mathbf{0},D)\,d\mathbf{b}_j
```

``` math
\ell(\theta)=\sum_{j=1}^G \log L_j(\theta)
```

### Laplace approximation used by `brsmm()`

Define
``` math
Q_j(\mathbf{b})=
\sum_{i=1}^{n_j}\log L_{ij}(\mathbf{b};\theta)+
\log\varphi_{q_b}(\mathbf{b};\mathbf{0},D)
```
and $`\hat{\mathbf{b}}_j=\arg\max_{\mathbf{b}} Q_j(\mathbf{b})`$, with
curvature
``` math
H_j = -\nabla^2 Q_j(\hat{\mathbf{b}}_j).
```
Then

``` math
\log L_j(\theta) \approx
Q_j(\hat{\mathbf{b}}_j) +
\frac{q_b}{2}\log(2\pi) -
\frac{1}{2}\log|H_j|.
```

[`brsmm()`](https://evandeilton.github.io/betaregscale/reference/brsmm.md)
maximizes the approximated $`\ell(\theta)`$ with
[`stats::optim()`](https://rdrr.io/r/stats/optim.html), and computes
group-level posterior modes $`\hat{\mathbf{b}}_j`$. For $`q_b = 1`$,
this reduces to the scalar random-intercept formula.

## Simulating clustered scale data

The next helper simulates data from a known mixed model to illustrate
fitting, inference, and recovery checks.

``` r

sim_brsmm_data <- function(seed = 3501L, g = 24L, ni = 12L,
                           beta = c(0.20, 0.65),
                           gamma = c(-0.15),
                           sigma_b = 0.55) {
  set.seed(seed)
  id <- factor(rep(seq_len(g), each = ni))
  n <- length(id)
  x1 <- rnorm(n)

  b_true <- rnorm(g, mean = 0, sd = sigma_b)
  eta_mu <- beta[1] + beta[2] * x1 + b_true[as.integer(id)]
  eta_phi <- rep(gamma[1], n)

  mu <- plogis(eta_mu)
  phi <- plogis(eta_phi)
  shp <- brs_repar(mu = mu, phi = phi, repar = 2)
  y <- round(stats::rbeta(n, shp$shape1, shp$shape2) * 100)

  list(
    data = data.frame(y = y, x1 = x1, id = id),
    truth = list(beta = beta, gamma = gamma, sigma_b = sigma_b, b = b_true)
  )
}

sim <- sim_brsmm_data(
  g = 12,
  ni = 20,
  beta = c(0.20, 0.65),
  gamma = c(-0.15),
  sigma_b = 0.55
)
str(sim$data)
#> 'data.frame':    240 obs. of  3 variables:
#>  $ y : num  0 83 99 56 65 10 8 1 65 98 ...
#>  $ x1: num  -0.3677 -2.0069 -0.0469 -0.2468 0.7634 ...
#>  $ id: Factor w/ 12 levels "1","2","3","4",..: 1 1 1 1 1 1 1 1 1 1 ...
```

## Fitting `brsmm()`

``` r

fit_mm <- brsmm(
  y ~ x1,
  random = ~ 1 | id,
  data = sim$data,
  repar = 2,
  int_method = "laplace",
  method = "BFGS",
  control = list(maxit = 1000)
)

summary(fit_mm)
#> 
#> Call:
#> brsmm(formula = y ~ x1, random = ~1 | id, data = sim$data, repar = 2, 
#>     int_method = "laplace", method = "BFGS", control = list(maxit = 1000))
#> 
#> Randomized Quantile Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -2.7630 -0.6319  0.0210  0.7206  3.5456 
#> 
#> Coefficients (mean model with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.36831    0.15411   2.390   0.0169 *  
#> x1           0.63299    0.09466   6.687 2.28e-11 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Phi coefficients (precision model with logit link):
#>             Estimate Std. Error z value Pr(>|z|)  
#> (Intercept) -0.15937    0.08481  -1.879   0.0602 .
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Random-effects parameters (Cholesky scale):
#>                      Estimate Std. Error z value Pr(>|z|)   
#> logSD.(Intercept)|id  -0.7974     0.2929  -2.723  0.00647 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Mixed beta interval model (Laplace)
#> Observations: 240  | Groups: 12 
#> Log-likelihood: -1008.2467 on 4 Df | AIC: 2024.4934 | BIC: 2038.4160 
#> Pseudo R-squared: 0.1364 
#> Number of iterations: 49 (BFGS) 
#> Censoring: 212 interval | 8 left | 20 right
```

## Random intercept + slope example

The model below includes a random intercept and random slope for `x1`:

``` r

fit_mm_rs <- brsmm(
  y ~ x1,
  random = ~ 1 + x1 | id,
  data = sim$data,
  repar = 2,
  int_method = "laplace",
  method = "BFGS",
  control = list(maxit = 1200)
)

summary(fit_mm_rs)
#> 
#> Call:
#> brsmm(formula = y ~ x1, random = ~1 + x1 | id, data = sim$data, 
#>     repar = 2, int_method = "laplace", method = "BFGS", control = list(maxit = 1200))
#> 
#> Randomized Quantile Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -2.9374 -0.6419  0.0104  0.6846  2.3799 
#> 
#> Coefficients (mean model with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)   0.3532     0.1579   2.238   0.0253 *  
#> x1            0.6292     0.1057   5.950 2.67e-09 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Phi coefficients (precision model with logit link):
#>             Estimate Std. Error z value Pr(>|z|)  
#> (Intercept) -0.17466    0.08499  -2.055   0.0399 *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Random-effects parameters (Cholesky scale):
#>                       Estimate Std. Error z value Pr(>|z|)   
#> logSD.(Intercept)|id   -0.7693     0.2867  -2.684  0.00728 **
#> cov.x1:(Intercept)|id  -0.1575     0.1079  -1.460  0.14421   
#> logSD.x1|id            -4.8447    10.8796  -0.445  0.65611   
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Mixed beta interval model (Laplace)
#> Observations: 240  | Groups: 12 
#> Log-likelihood: -1007.0654 on 6 Df | AIC: 2026.1308 | BIC: 2047.0146 
#> Pseudo R-squared: 0.1364 
#> Number of iterations: 32 (BFGS) 
#> Censoring: 212 interval | 8 left | 20 right
```

Covariance structure of random effects:

``` r

kbl10(fit_mm_rs$random$D)
```

|   V1    |   V2    |
|:-------:|:-------:|
| 0.2147  | -0.0730 |
| -0.0730 | 0.0249  |

``` r

kbl10(
  data.frame(term = names(fit_mm_rs$random$sd_b), sd = as.numeric(fit_mm_rs$random$sd_b)),
  digits = 4
)
```

|    term     |   sd   |
|:-----------:|:------:|
| (Intercept) | 0.4633 |
|     x1      | 0.1577 |

``` r

kbl10(head(ranef(fit_mm_rs), 10))
```

| (Intercept) |   x1    |
|:-----------:|:-------:|
|   -0.0039   | 0.0012  |
|   0.2585    | -0.0880 |
|   0.0431    | -0.0146 |
|   0.1528    | -0.0518 |
|   -0.0353   | 0.0119  |
|   0.6181    | -0.2101 |
|   -0.4099   | 0.1396  |
|   -0.1081   | 0.0365  |
|   -0.4718   | 0.1605  |
|   -0.6701   | 0.2278  |

## Additional studies of random effects (numerical and visual)

Following practices from established mixed-models packages, the package
now allows for a dedicated study of the random effects focusing on:

- $`D`$ structure and correlation;
- empirical distribution of modes by group;
- empirical shrinkage intensity;
- specific visual diagnostics for the random components.

``` r

re_study <- brsmm_re_study(fit_mm_rs)
print(re_study)
#> 
#> Random-effects study
#> Groups: 12 
#> 
#> Random-effects (VarCorr):
#>   Name                      Std.Dev.  Corr
#>   re1                         0.4633
#>   re2                         0.1577  -0.9988
#> 
#> ICC (latent logistic scale): 0.0613
#> 
#> Summary by term (SD_model = model SD; shrinkage = Var(modes)/Var(model)):
#>         term sd_model mean_mode sd_mode shrinkage_ratio shapiro_p
#>  (Intercept)   0.4633     6e-04  0.4159          0.8056    0.8381
#>           x1   0.1577    -2e-04  0.1414          0.8039    0.8370
kbl10(re_study$summary)
```

|    term     | sd_model | mean_mode | sd_mode | shrinkage_ratio | shapiro_p |
|:-----------:|:--------:|:---------:|:-------:|:---------------:|:---------:|
| (Intercept) |  0.4633  |   6e-04   | 0.4159  |     0.8056      |  0.8381   |
|     x1      |  0.1577  |  -2e-04   | 0.1414  |     0.8039      |  0.8370   |

``` r

kbl10(re_study$D)
```

|   V1    |   V2    |
|:-------:|:-------:|
| 0.2147  | -0.0730 |
| -0.0730 | 0.0249  |

``` r

kbl10(re_study$Corr)
```

|   V1    |   V2    |
|:-------:|:-------:|
| 1.0000  | -0.9988 |
| -0.9988 | 1.0000  |

Suggested visualizations for random effects:

``` r

if (requireNamespace("ggplot2", quietly = TRUE)) {
  autoplot.brsmm(fit_mm_rs, type = "ranef_caterpillar")
  autoplot.brsmm(fit_mm_rs, type = "ranef_density")
  autoplot.brsmm(fit_mm_rs, type = "ranef_pairs")
  autoplot.brsmm(fit_mm_rs, type = "ranef_qq")
}
```

![](brs-mm_files/figure-html/ranef-visuals-1.png)

## Core methods

### Coefficients and random effects

`coef(fit_mm, model = "random")` returns packed random-effect covariance
parameters on the optimizer scale (lower-Cholesky, with a log-diagonal).
For random-intercept models, this simplifies to $`\log \sigma_b`$.

``` r

kbl10(
  data.frame(
    parameter = names(coef(fit_mm, model = "full")),
    estimate = as.numeric(coef(fit_mm, model = "full"))
  ),
  digits = 4
)
```

|            parameter             | estimate |
|:--------------------------------:|:--------:|
|           (Intercept)            |  0.3683  |
|                x1                |  0.6330  |
|        (phi)\_(Intercept)        | -0.1594  |
| (re_chol_logsd)\_(Intercept)\|id | -0.7974  |

``` r

kbl10(
  data.frame(
    log_sigma_b = as.numeric(coef(fit_mm, model = "random")),
    sigma_b = as.numeric(exp(coef(fit_mm, model = "random")))
  ),
  digits = 4
)
```

| log_sigma_b | sigma_b |
|:-----------:|:-------:|
|   -0.7974   | 0.4505  |

``` r

kbl10(head(ranef(fit_mm), 10))
```

|    x    |
|:-------:|
| -0.0420 |
| 0.2152  |
| 0.0468  |
| 0.2022  |
| -0.0720 |
| 0.6097  |
| -0.3624 |
| -0.1838 |
| -0.4129 |
| -0.6104 |

For random intercept + slope models:

``` r

kbl10(
  data.frame(
    parameter = names(coef(fit_mm_rs, model = "random")),
    estimate = as.numeric(coef(fit_mm_rs, model = "random"))
  ),
  digits = 4
)
```

|            parameter             | estimate |
|:--------------------------------:|:--------:|
| (re_chol_logsd)\_(Intercept)\|id | -0.7693  |
|  (re_chol)\_x1:(Intercept)\|id   | -0.1575  |
|     (re_chol_logsd)\_x1\|id      | -4.8447  |

``` r

kbl10(fit_mm_rs$random$D)
```

|   V1    |   V2    |
|:-------:|:-------:|
| 0.2147  | -0.0730 |
| -0.0730 | 0.0249  |

### Variance-covariance, summary and likelihood criteria

``` r

vc <- vcov(fit_mm)
dim(vc)
#> [1] 4 4

sm <- summary(fit_mm)
kbl10(sm$coefficients)
```

|  | mean.Estimate | mean.Std..Error | mean.z.value | mean.Pr…z.. | precision.Estimate | precision.Std..Error | precision.z.value | precision.Pr…z.. | random.Estimate | random.Std..Error | random.z.value | random.Pr…z.. |
|:---|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| (Intercept) | 0.3683 | 0.1541 | 2.3898 | 0.0169 | -0.1594 | 0.0848 | -1.8791 | 0.0602 | -0.7974 | 0.2929 | -2.7228 | 0.0065 |
| x1 | 0.6330 | 0.0947 | 6.6869 | 0.0000 | -0.1594 | 0.0848 | -1.8791 | 0.0602 | -0.7974 | 0.2929 | -2.7228 | 0.0065 |

``` r


kbl10(
  data.frame(
    logLik = as.numeric(logLik(fit_mm)),
    AIC = AIC(fit_mm),
    BIC = BIC(fit_mm),
    nobs = nobs(fit_mm)
  ),
  digits = 4
)
```

|  logLik   |   AIC    |   BIC    | nobs |
|:---------:|:--------:|:--------:|:----:|
| -1008.247 | 2024.493 | 2038.416 | 240  |

### Fitted values, prediction and residuals

``` r

kbl10(
  data.frame(
    mu_hat = head(fitted(fit_mm, type = "mu")),
    phi_hat = head(fitted(fit_mm, type = "phi")),
    pred_mu = head(predict(fit_mm, type = "response")),
    pred_eta = head(predict(fit_mm, type = "link")),
    pred_phi = head(predict(fit_mm, type = "precision")),
    pred_var = head(predict(fit_mm, type = "variance"))
  ),
  digits = 4
)
```

| mu_hat | phi_hat | pred_mu | pred_eta | pred_phi | pred_var |
|:------:|:-------:|:-------:|:--------:|:--------:|:--------:|
| 0.5234 | 0.4602  | 0.5234  |  0.0936  |  0.4602  |  0.1148  |
| 0.2801 | 0.4602  | 0.2801  | -0.9440  |  0.4602  |  0.0928  |
| 0.5736 | 0.4602  | 0.5736  |  0.2966  |  0.4602  |  0.1126  |
| 0.5424 | 0.4602  | 0.5424  |  0.1701  |  0.4602  |  0.1142  |
| 0.6920 | 0.4602  | 0.6920  |  0.8096  |  0.4602  |  0.0981  |
| 0.4792 | 0.4602  | 0.4792  | -0.0834  |  0.4602  |  0.1149  |

``` r


kbl10(
  data.frame(
    res_response = head(residuals(fit_mm, type = "response")),
    res_pearson = head(residuals(fit_mm, type = "pearson"))
  ),
  digits = 4
)
```

| res_response | res_pearson |
|:------------:|:-----------:|
|   -0.5234    |   -1.5446   |
|    0.5499    |   1.8052    |
|    0.4164    |   1.2410    |
|    0.0176    |   0.0520    |
|   -0.0420    |   -0.1342   |
|   -0.3792    |   -1.1188   |

### Diagnostic plotting methods

[`plot.brsmm()`](https://evandeilton.github.io/betaregscale/reference/plot.brsmm.md)
supports base and ggplot2 backends:

``` r

plot(fit_mm, which = 1:4, type = "pearson")
```

![](brs-mm_files/figure-html/methods-plot-1.png)

``` r


if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot(fit_mm, which = 1:2, gg = TRUE)
}
```

![](brs-mm_files/figure-html/methods-plot-2.png)

[`autoplot.brsmm()`](https://evandeilton.github.io/betaregscale/reference/autoplot.brsmm.md)
provides focused ggplot diagnostics:

``` r

if (requireNamespace("ggplot2", quietly = TRUE)) {
  autoplot.brsmm(fit_mm, type = "calibration")
  autoplot.brsmm(fit_mm, type = "score_dist")
  autoplot.brsmm(fit_mm, type = "ranef_qq")
  autoplot.brsmm(fit_mm, type = "residuals_by_group")
}
```

![](brs-mm_files/figure-html/methods-autoplot-1.png)

### Prediction with `newdata`

If `newdata` contains unseen groups,
[`predict.brsmm()`](https://evandeilton.github.io/betaregscale/reference/predict.brsmm.md)
uses a random effect equal to zero for those levels.

``` r

nd <- sim$data[1:8, c("x1", "id")]
kbl10(
  data.frame(pred_seen = as.numeric(predict(fit_mm, newdata = nd, type = "response"))),
  digits = 4
)
```

| pred_seen |
|:---------:|
|  0.5234   |
|  0.2801   |
|  0.5736   |
|  0.5424   |
|  0.6920   |
|  0.4792   |
|  0.6047   |
|  0.4668   |

``` r


nd_unseen <- nd
nd_unseen$id <- factor(rep("new_cluster", nrow(nd_unseen)))
kbl10(
  data.frame(pred_unseen = as.numeric(predict(fit_mm, newdata = nd_unseen, type = "response"))),
  digits = 4
)
```

| pred_unseen |
|:-----------:|
|   0.5338    |
|   0.2886    |
|   0.5839    |
|   0.5528    |
|   0.7009    |
|   0.4897    |
|   0.6147    |
|   0.4773    |

The same logic applies to random intercept + slope models:

``` r

kbl10(
  data.frame(pred_rs_seen = as.numeric(predict(fit_mm_rs, newdata = nd, type = "response"))),
  digits = 4
)
```

| pred_rs_seen |
|:------------:|
|    0.5293    |
|    0.2858    |
|    0.5793    |
|    0.5483    |
|    0.6965    |
|    0.4853    |
|    0.6101    |
|    0.4730    |

``` r

kbl10(
  data.frame(pred_rs_unseen = as.numeric(predict(fit_mm_rs, newdata = nd_unseen, type = "response"))),
  digits = 4
)
```

| pred_rs_unseen |
|:--------------:|
|     0.5304     |
|     0.2871     |
|     0.5802     |
|     0.5493     |
|     0.6971     |
|     0.4865     |
|     0.6109     |
|     0.4742     |

## Statistical tests and validation workflow

### Wald tests (from `summary`)

[`summary.brsmm()`](https://evandeilton.github.io/betaregscale/reference/summary.brsmm.md)
reports Wald $`z`$-tests for each parameter:
``` math
z_k = \hat\theta_k / \mathrm{SE}(\hat\theta_k).
```

``` r

sm <- summary(fit_mm)
kbl10(sm$coefficients)
```

|  | mean.Estimate | mean.Std..Error | mean.z.value | mean.Pr…z.. | precision.Estimate | precision.Std..Error | precision.z.value | precision.Pr…z.. | random.Estimate | random.Std..Error | random.z.value | random.Pr…z.. |
|:---|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| (Intercept) | 0.3683 | 0.1541 | 2.3898 | 0.0169 | -0.1594 | 0.0848 | -1.8791 | 0.0602 | -0.7974 | 0.2929 | -2.7228 | 0.0065 |
| x1 | 0.6330 | 0.0947 | 6.6869 | 0.0000 | -0.1594 | 0.0848 | -1.8791 | 0.0602 | -0.7974 | 0.2929 | -2.7228 | 0.0065 |

### Evolutionary scheme and Likelihood Ratio (LR) test selection

A practical workflow of increasing complexity:

1.  [`brs()`](https://evandeilton.github.io/betaregscale/reference/brs.md):
    no random effect (ignores clustering);
2.  `brsmm(..., random = ~ 1 | id)`: random intercept;
3.  `brsmm(..., random = ~ 1 + x1 | id)`: random intercept + slope.

In the first jump (`brs` to `brsmm` with intercept), the hypothesis
$`\sigma_b^2 = 0`$ lies on the boundary of the parameter space. Thus,
the classical asymptotic $`\chi^2`$ reference distribution should be
interpreted with caution. In the second jump (intercept to intercept +
slope), the Likelihood Ratio (LR) test with a $`\chi^2`$ distribution is
commonly used as a practical diagnostic for goodness-of-fit gains.

``` r

# Base model without a random effect
fit_brs <- brs(
  y ~ x1,
  data = sim$data,
  repar = 2
)

# Reuse the mixed models already fitted:
# fit_mm    : random = ~ 1 | id
# fit_mm_rs : random = ~ 1 + x1 | id

tab_lr <- anova(fit_brs, fit_mm, fit_mm_rs, test = "Chisq")
kbl10(
  data.frame(model = rownames(tab_lr), tab_lr, row.names = NULL),
  digits = 4
)
```

|   model    | Df  |  logLik   |   AIC    |   BIC    |  Chisq  | Chi.Df | Pr..Chisq. |
|:----------:|:---:|:---------:|:--------:|:--------:|:-------:|:------:|:----------:|
|  M1 (brs)  |  3  | -1014.913 | 2035.826 | 2046.268 |   NA    |   NA   |     NA     |
| M2 (brsmm) |  4  | -1008.247 | 2024.493 | 2038.416 | 13.3331 |   1    |   0.0003   |
| M3 (brsmm) |  6  | -1007.065 | 2026.131 | 2047.015 | 2.3627  |   2    |   0.3069   |

Operational decision rule (analytical):

- If the second jump (RI to RI+RS) does not improve the fit (high
  p-value), prefer the random-intercept model for parsimony.
- If there is a robust gain, adopt the RI+RS model and validate
  parameter stability (especially `sd_b` and the $`D`$ matrix) via
  sensitivity and residual diagnostics.

### Residual diagnostics (quick checks)

``` r

r <- residuals(fit_mm, type = "pearson")
kbl10(
  data.frame(
    mean = mean(r),
    sd = stats::sd(r),
    q025 = as.numeric(stats::quantile(r, 0.025)),
    q975 = as.numeric(stats::quantile(r, 0.975))
  ),
  digits = 4
)
```

|  mean  |   sd   |  q025  |  q975  |
|:------:|:------:|:------:|:------:|
| 0.0305 | 0.9788 | -1.842 | 1.5132 |

## Parameter recovery experiment

A single-fit recovery table can be produced directly from the previous
fit:

``` r

est <- c(
  beta0 = unname(coef(fit_mm, model = "mean")[1]),
  beta1 = unname(coef(fit_mm, model = "mean")[2]),
  sigma_b = unname(exp(coef(fit_mm, model = "random")))
)

true <- c(
  beta0 = sim$truth$beta[1],
  beta1 = sim$truth$beta[2],
  sigma_b = sim$truth$sigma_b
)

recovery_table <- data.frame(
  parameter = names(true),
  true = as.numeric(true),
  estimate = as.numeric(est[names(true)]),
  bias = as.numeric(est[names(true)] - true)
)
kbl10(recovery_table)
```

| parameter | true | estimate |  bias   |
|:---------:|:----:|:--------:|:-------:|
|   beta0   | 0.20 |  0.3683  | 0.1683  |
|   beta1   | 0.65 |  0.6330  | -0.0170 |
|  sigma_b  | 0.55 |  0.4505  | -0.0995 |

For a Monte Carlo recovery study, repeat simulation and fitting across
replicates:

``` r

mc_recovery <- function(R = 50L, seed = 7001L) {
  set.seed(seed)
  out <- vector("list", R)

  for (r in seq_len(R)) {
    sim_r <- sim_brsmm_data(seed = seed + r)
    fit_r <- brsmm(
      y ~ x1,
      random = ~ 1 | id,
      data = sim_r$data,
      repar = 2,
      int_method = "laplace",
      method = "BFGS",
      control = list(maxit = 1000)
    )

    out[[r]] <- c(
      beta0 = unname(coef(fit_r, model = "mean")[1]),
      beta1 = unname(coef(fit_r, model = "mean")[2]),
      sigma_b = unname(exp(coef(fit_r, model = "random")))
    )
  }

  est <- do.call(rbind, out)
  truth <- c(beta0 = 0.20, beta1 = 0.65, sigma_b = 0.55)

  data.frame(
    parameter = colnames(est),
    truth = as.numeric(truth[colnames(est)]),
    mean_est = colMeans(est),
    bias = colMeans(est) - truth[colnames(est)],
    rmse = sqrt(colMeans((sweep(est, 2, truth[colnames(est)], "-"))^2))
  )
}

kbl10(mc_recovery(R = 50))
```

## How this maps to automated package tests

The package test suite includes dedicated `brsmm` tests for:

1.  fitting with Laplace integration;
2.  one- and two-part formulas;
3.  S3 methods (`coef`, `vcov`, `summary`, `predict`, `residuals`,
    `ranef`);
4.  parameter recovery under known DGP settings.

Run locally:

``` r

devtools::test(filter = "brsmm")
```

## References

Ferrari, S. L. P. and Cribari-Neto, F. (2004). Beta regression for
modelling rates and proportions. *Journal of Applied Statistics*, 31(7),
799-815. DOI: 10.1080/0266476042000214501. Validated online via:
<https://doi.org/10.1080/0266476042000214501>.

Pinheiro, J. C. and Bates, D. M. (2000). *Mixed-Effects Models in S and
S-PLUS*. Springer. DOI: 10.1007/b98882. Validated online via:
<https://doi.org/10.1007/b98882>.

Rue, H., Martino, S., and Chopin, N. (2009). Approximate Bayesian
inference for latent Gaussian models by using integrated nested Laplace
approximations. *Journal of the Royal Statistical Society: Series B*,
71(2), 319-392. DOI: 10.1111/j.1467-9868.2008.00700.x. Validated online
via: <https://doi.org/10.1111/j.1467-9868.2008.00700.x>.
