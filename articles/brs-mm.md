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

Assume observations $i = 1,\ldots,n_{j}$ within groups $j = 1,\ldots,G$,
with group-specific random-effects vector
$\mathbf{b}_{j} \in {\mathbb{R}}^{q_{b}}$.

### Linear predictors

$$\eta_{\mu,ij} = x_{ij}^{\top}\beta + w_{ij}^{\top}\mathbf{b}_{j},\qquad\eta_{\phi,ij} = z_{ij}^{\top}\gamma$$

$$\mu_{ij} = g^{- 1}\left( \eta_{\mu,ij} \right),\qquad\phi_{ij} = h^{- 1}\left( \eta_{\phi,ij} \right)$$

with $g( \cdot )$ and $h( \cdot )$ chosen by `link` and `link_phi`. The
random-effects design row $w_{ij}$ is defined by
`random = ~ terms | group`.

### Beta parameterization

For each $\left( \mu_{ij},\phi_{ij} \right)$, `repar` maps to beta shape
parameters $\left( a_{ij},b_{ij} \right)$ via
[`brs_repar()`](https://evandeilton.github.io/betaregscale/reference/brs_repar.md).

### Conditional contribution by censoring type

Each observation contributes:

$$L_{ij}\left( b_{j};\theta \right) = \begin{cases}
{f\left( y_{ij};a_{ij},b_{ij} \right),} & {\delta_{ij} = 0} \\
{F\left( u_{ij};a_{ij},b_{ij} \right),} & {\delta_{ij} = 1} \\
{1 - F\left( l_{ij};a_{ij},b_{ij} \right),} & {\delta_{ij} = 2} \\
{F\left( u_{ij};a_{ij},b_{ij} \right) - F\left( l_{ij};a_{ij},b_{ij} \right),} & {\delta_{ij} = 3}
\end{cases}$$

where $l_{ij},u_{ij}$ are interval endpoints on $(0,1)$, $f( \cdot )$ is
beta density, and $F( \cdot )$ is beta CDF.

### Random-effects distribution

$$\mathbf{b}_{j} \sim \mathcal{N}(\mathbf{0},D),$$

where $D$ is a symmetric positive-definite covariance matrix.
Internally,
[`brsmm()`](https://evandeilton.github.io/betaregscale/reference/brsmm.md)
optimizes a packed lower-Cholesky parameterization $D = LL^{\top}$
(diagonal entries on log-scale for positivity).

### Group marginal likelihood

$$L_{j}(\theta) = \int_{{\mathbb{R}}^{q_{b}}}\prod\limits_{i = 1}^{n_{j}}L_{ij}\left( b_{j};\theta \right)\,\varphi_{q_{b}}\left( \mathbf{b}_{j};\mathbf{0},D \right)\, d\mathbf{b}_{j}$$

$$\ell(\theta) = \sum\limits_{j = 1}^{G}\log L_{j}(\theta)$$

### Laplace approximation used by `brsmm()`

Define
$$Q_{j}(\mathbf{b}) = \sum\limits_{i = 1}^{n_{j}}\log L_{ij}(\mathbf{b};\theta) + \log\varphi_{q_{b}}(\mathbf{b};\mathbf{0},D)$$
and
${\widehat{\mathbf{b}}}_{j} = \arg\max_{\mathbf{b}}Q_{j}(\mathbf{b})$,
with curvature
$$H_{j} = - \nabla^{2}Q_{j}\left( {\widehat{\mathbf{b}}}_{j} \right).$$
Then

$$\log L_{j}(\theta) \approx Q_{j}\left( {\widehat{\mathbf{b}}}_{j} \right) + \frac{q_{b}}{2}\log(2\pi) - \frac{1}{2}\log\left| H_{j} \right|.$$

[`brsmm()`](https://evandeilton.github.io/betaregscale/reference/brsmm.md)
maximizes the approximated $\ell(\theta)$ with
[`stats::optim()`](https://rdrr.io/r/stats/optim.html), and computes
group-level posterior modes ${\widehat{\mathbf{b}}}_{j}$. For
$q_{b} = 1$, this reduces to the scalar random-intercept formula.

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

sim <- sim_brsmm_data()
str(sim$data)
#> 'data.frame':    288 obs. of  3 variables:
#>  $ y : num  18 29 12 43 56 55 3 21 4 4 ...
#>  $ x1: num  -0.3677 -2.0069 -0.0469 -0.2468 0.7634 ...
#>  $ id: Factor w/ 24 levels "1","2","3","4",..: 1 1 1 1 1 1 1 1 1 1 ...
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
#> Mixed beta interval model (Laplace)
#> Observations: 288  | Groups: 24 
#> logLik =-1231.3848 | AIC =2470.7696 | BIC =2485.4215
#> 
#>                                Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)                     0.01050    0.07457   0.141  0.88802    
#> x1                              0.62964    0.08652   7.278  3.4e-13 ***
#> (phi)_(Intercept)              -0.13629    0.07738  -1.761  0.07819 .  
#> (re_chol_logsd)_(Intercept)|id -0.84336    0.26108  -3.230  0.00124 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

## Random intercept + slope example

The model below includes random intercept and random slope for `x1`:

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
#> Mixed beta interval model (Laplace)
#> Observations: 288  | Groups: 24 
#> logLik =-1229.9565 | AIC =2471.9130 | BIC =2493.8907
#> 
#>                                 Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)                     0.005468   0.040347   0.136  0.89220    
#> x1                              0.633252   0.107380   5.897 3.69e-09 ***
#> (phi)_(Intercept)              -0.159369   0.079329  -2.009  0.04454 *  
#> (re_chol_logsd)_(Intercept)|id -0.904046   0.283596  -3.188  0.00143 ** 
#> (re_chol)_x1:(Intercept)|id     0.031810   0.088554   0.359  0.71943    
#> (re_chol_logsd)_x1|id          -1.248798   0.423772  -2.947  0.00321 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Covariance structure of random effects:

``` r
knitr::kable(fit_mm_rs$random$D, digits = 4)
```

|        |        |
|-------:|-------:|
| 0.1640 | 0.0129 |
| 0.0129 | 0.0833 |

``` r
knitr::kable(
  data.frame(term = names(fit_mm_rs$random$sd_b), sd = as.numeric(fit_mm_rs$random$sd_b)),
  digits = 4
)
```

| term        |     sd |
|:------------|-------:|
| (Intercept) | 0.4049 |
| x1          | 0.2886 |

``` r
knitr::kable(head(ranef.brsmm(fit_mm_rs)), digits = 4)
```

| (Intercept) |      x1 |
|------------:|--------:|
|     -0.3977 | -0.2520 |
|      0.4303 | -0.0387 |
|     -0.1598 |  0.0259 |
|     -0.3099 |  0.3768 |
|      0.0120 | -0.1438 |
|      0.3996 | -0.1521 |

## Estudos adicionais dos efeitos aleatórios (numéricos e visuais)

Seguindo práticas de pacotes mistos consolidados, o pacote agora permite
um estudo dedicado dos efeitos aleatórios com foco em:

- estrutura $D$ e correlação;
- distribuição empírica dos modos por grupo;
- intensidade de shrinkage empírico;
- diagnósticos visuais específicos para os componentes aleatórios.

``` r
re_study <- brsmm_re_study(fit_mm_rs)
print(re_study)
#> 
#> Random-effects study
#> Groups: 24 
#> 
#> Summary by term:
#>         term sd_model mean_mode sd_mode shrinkage_ratio shapiro_p
#>  (Intercept)   0.4049   -0.0012  0.3044          0.5653    0.0543
#>           x1   0.2886   -0.0081  0.1686          0.3412    0.0199
#> 
#> Estimated covariance matrix D:
#>        [,1]   [,2]
#> [1,] 0.1640 0.0129
#> [2,] 0.0129 0.0833
#> 
#> Estimated correlation matrix:
#>        [,1]   [,2]
#> [1,] 1.0000 0.1102
#> [2,] 0.1102 1.0000
knitr::kable(re_study$summary, digits = 4)
```

| term        | sd_model | mean_mode | sd_mode | shrinkage_ratio | shapiro_p |
|:------------|---------:|----------:|--------:|----------------:|----------:|
| (Intercept) |   0.4049 |   -0.0012 |  0.3044 |          0.5653 |    0.0543 |
| x1          |   0.2886 |   -0.0081 |  0.1686 |          0.3412 |    0.0199 |

``` r
knitr::kable(re_study$D, digits = 4)
```

|        |        |
|-------:|-------:|
| 0.1640 | 0.0129 |
| 0.0129 | 0.0833 |

``` r
knitr::kable(re_study$Corr, digits = 4)
```

|        |        |
|-------:|-------:|
| 1.0000 | 0.1102 |
| 0.1102 | 1.0000 |

Visualizações sugeridas para efeitos aleatórios:

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
parameters on optimizer scale (lower-Cholesky, with log-diagonal). For
random-intercept models this simplifies to $\log\sigma_{b}$.

``` r
knitr::kable(
  data.frame(parameter = names(coef(fit_mm, model = "full")),
             estimate = as.numeric(coef(fit_mm, model = "full"))),
  digits = 4
)
```

| parameter                        | estimate |
|:---------------------------------|---------:|
| (Intercept)                      |   0.0105 |
| x1                               |   0.6296 |
| (phi)\_(Intercept)               |  -0.1363 |
| (re_chol_logsd)\_(Intercept)\|id |  -0.8434 |

``` r
knitr::kable(
  data.frame(log_sigma_b = as.numeric(coef(fit_mm, model = "random")),
             sigma_b = as.numeric(exp(coef(fit_mm, model = "random")))),
  digits = 4
)
```

| log_sigma_b | sigma_b |
|------------:|--------:|
|     -0.8434 |  0.4303 |

``` r
knitr::kable(head(ranef.brsmm(fit_mm)), digits = 4)
```

|       x |
|--------:|
| -0.3708 |
|  0.4555 |
| -0.1742 |
| -0.4594 |
|  0.0248 |
|  0.4296 |

For random intercept + slope models:

``` r
knitr::kable(
  data.frame(
    parameter = names(coef(fit_mm_rs, model = "random")),
    estimate = as.numeric(coef(fit_mm_rs, model = "random"))
  ),
  digits = 4
)
```

| parameter                        | estimate |
|:---------------------------------|---------:|
| (re_chol_logsd)\_(Intercept)\|id |  -0.9040 |
| (re_chol)\_x1:(Intercept)\|id    |   0.0318 |
| (re_chol_logsd)\_x1\|id          |  -1.2488 |

``` r
knitr::kable(fit_mm_rs$random$D, digits = 4)
```

|        |        |
|-------:|-------:|
| 0.1640 | 0.0129 |
| 0.0129 | 0.0833 |

### Variance-covariance, summary and likelihood criteria

``` r
vc <- vcov(fit_mm)
dim(vc)
#> [1] 4 4

sm <- summary(fit_mm)
knitr::kable(sm$coefficients, digits = 4)
```

|                                  | Estimate | Std. Error | z value | Pr(\>\|z\|) |
|:---------------------------------|---------:|-----------:|--------:|------------:|
| (Intercept)                      |   0.0105 |     0.0746 |  0.1408 |      0.8880 |
| x1                               |   0.6296 |     0.0865 |  7.2777 |      0.0000 |
| (phi)\_(Intercept)               |  -0.1363 |     0.0774 | -1.7613 |      0.0782 |
| (re_chol_logsd)\_(Intercept)\|id |  -0.8434 |     0.2611 | -3.2303 |      0.0012 |

``` r

knitr::kable(
  data.frame(
    logLik = as.numeric(logLik(fit_mm)),
    AIC = AIC(fit_mm),
    BIC = BIC(fit_mm),
    nobs = nobs(fit_mm)
  ),
  digits = 4
)
```

|    logLik |     AIC |      BIC | nobs |
|----------:|--------:|---------:|-----:|
| -1231.385 | 2470.77 | 2485.421 |  288 |

### Fitted values, prediction and residuals

``` r
knitr::kable(
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
|-------:|--------:|--------:|---------:|---------:|---------:|
| 0.3562 |   0.466 |  0.3562 |  -0.5918 |    0.466 |   0.1069 |
| 0.1647 |   0.466 |  0.1647 |  -1.6239 |    0.466 |   0.0641 |
| 0.4038 |   0.466 |  0.4038 |  -0.3898 |    0.466 |   0.1122 |
| 0.3739 |   0.466 |  0.3739 |  -0.5157 |    0.466 |   0.1091 |
| 0.5301 |   0.466 |  0.5301 |   0.1204 |    0.466 |   0.1161 |
| 0.3170 |   0.466 |  0.3170 |  -0.7678 |    0.466 |   0.1009 |

``` r

knitr::kable(
  data.frame(
    res_response = head(residuals(fit_mm, type = "response")),
    res_pearson = head(residuals(fit_mm, type = "pearson"))
  ),
  digits = 4
)
```

| res_response | res_pearson |
|-------------:|------------:|
|      -0.1762 |     -0.5391 |
|       0.1253 |      0.4950 |
|      -0.2838 |     -0.8472 |
|       0.0561 |      0.1700 |
|       0.0299 |      0.0879 |
|       0.2330 |      0.7337 |

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
uses random effect equal to zero for those levels.

``` r
nd <- sim$data[1:8, c("x1", "id")]
knitr::kable(
  data.frame(pred_seen = as.numeric(predict(fit_mm, newdata = nd, type = "response"))),
  digits = 4
)
```

| pred_seen |
|----------:|
|    0.3562 |
|    0.1647 |
|    0.4038 |
|    0.3739 |
|    0.5301 |
|    0.3170 |
|    0.4348 |
|    0.3064 |

``` r

nd_unseen <- nd
nd_unseen$id <- factor(rep("new_cluster", nrow(nd_unseen)))
knitr::kable(
  data.frame(pred_unseen = as.numeric(predict(fit_mm, newdata = nd_unseen, type = "response"))),
  digits = 4
)
```

| pred_unseen |
|------------:|
|      0.4450 |
|      0.2222 |
|      0.4952 |
|      0.4638 |
|      0.6204 |
|      0.4020 |
|      0.5271 |
|      0.3902 |

The same logic applies to random intercept + slope models:

``` r
knitr::kable(
  data.frame(pred_rs_seen = as.numeric(predict(fit_mm_rs, newdata = nd, type = "response"))),
  digits = 4
)
```

| pred_rs_seen |
|-------------:|
|       0.3699 |
|       0.2391 |
|       0.3989 |
|       0.3808 |
|       0.4747 |
|       0.3455 |
|       0.4175 |
|       0.3388 |

``` r
knitr::kable(
  data.frame(pred_rs_unseen = as.numeric(predict(fit_mm_rs, newdata = nd_unseen, type = "response"))),
  digits = 4
)
```

| pred_rs_unseen |
|---------------:|
|         0.4434 |
|         0.2200 |
|         0.4939 |
|         0.4624 |
|         0.6198 |
|         0.4003 |
|         0.5260 |
|         0.3884 |

## Statistical tests and validation workflow

### Wald tests (from `summary`)

[`summary.brsmm()`](https://evandeilton.github.io/betaregscale/reference/summary.brsmm.md)
reports Wald $z$-tests for each parameter:
$$z_{k} = {\widehat{\theta}}_{k}/{SE}\left( {\widehat{\theta}}_{k} \right).$$

``` r
sm <- summary(fit_mm)
knitr::kable(sm$coefficients, digits = 4)
```

|                                  | Estimate | Std. Error | z value | Pr(\>\|z\|) |
|:---------------------------------|---------:|-----------:|--------:|------------:|
| (Intercept)                      |   0.0105 |     0.0746 |  0.1408 |      0.8880 |
| x1                               |   0.6296 |     0.0865 |  7.2777 |      0.0000 |
| (phi)\_(Intercept)               |  -0.1363 |     0.0774 | -1.7613 |      0.0782 |
| (re_chol_logsd)\_(Intercept)\|id |  -0.8434 |     0.2611 | -3.2303 |      0.0012 |

### Esquema evolutivo e escolha por teste LR

Um fluxo prático de complexidade crescente:

1.  [`brs()`](https://evandeilton.github.io/betaregscale/reference/brs.md):
    sem efeito aleatório (ignora agrupamento);
2.  `brsmm(..., random = ~ 1 | id)`: intercepto aleatório;
3.  `brsmm(..., random = ~ 1 + x1 | id)`: intercepto + inclinação
    aleatórios.

No primeiro salto (`brs` $\rightarrow$`brsmm` com intercepto), a
hipótese $\sigma_{b}^{2} = 0$ está na fronteira do espaço paramétrico.
Assim, o referencial assintótico clássico $\chi^{2}$ deve ser
interpretado com cautela. No segundo salto (intercepto $\rightarrow$
intercepto + inclinação), o LR com $\chi^{2}$ costuma ser usado como
diagnóstico prático de ganho de ajuste.

``` r
# Modelo base sem efeito aleatório
fit_brs <- brs(
  y ~ x1,
  data = sim$data,
  repar = 2
)

# Reutiliza os ajustes mistos já estimados:
# fit_mm    : random = ~ 1 | id
# fit_mm_rs : random = ~ 1 + x1 | id

tab_lr <- anova(fit_brs, fit_mm, fit_mm_rs, test = "Chisq")
knitr::kable(
  data.frame(model = rownames(tab_lr), tab_lr, row.names = NULL),
  digits = 4
)
```

| model      |  Df |    logLik |      AIC |      BIC |   Chisq | Chi.Df | Pr..Chisq. |
|:-----------|----:|----------:|---------:|---------:|--------:|-------:|-----------:|
| M1 (brs)   |   3 | -1236.836 | 2479.672 | 2490.661 |      NA |     NA |         NA |
| M2 (brsmm) |   4 | -1231.385 | 2470.770 | 2485.421 | 10.9028 |      1 |     0.0010 |
| M3 (brsmm) |   6 | -1229.957 | 2471.913 | 2493.891 |  2.8567 |      2 |     0.2397 |

Regra operacional de decisão (analítica):

- Se o segundo salto (`RI -> RI+RS`) não melhora o ajuste (p alto),
  prefira o modelo com intercepto aleatório por parcimônia.
- Se houver ganho robusto, adote `RI+RS` e valide estabilidade dos
  parâmetros (especialmente `sd_b` e matriz `D`) por sensibilidade e
  diagnóstico residual.

### Residual diagnostics (quick checks)

``` r
r <- residuals(fit_mm, type = "pearson")
knitr::kable(
  data.frame(
    mean = mean(r),
    sd = stats::sd(r),
    q025 = as.numeric(stats::quantile(r, 0.025)),
    q975 = as.numeric(stats::quantile(r, 0.975))
  ),
  digits = 4
)
```

|  mean |     sd |    q025 |   q975 |
|------:|-------:|--------:|-------:|
| 0.038 | 0.9611 | -1.5667 | 1.6286 |

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
knitr::kable(recovery_table, digits = 4)
```

| parameter | true | estimate |    bias |
|:----------|-----:|---------:|--------:|
| beta0     | 0.20 |   0.0105 | -0.1895 |
| beta1     | 0.65 |   0.6296 | -0.0204 |
| sigma_b   | 0.55 |   0.4303 | -0.1197 |

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

mc_recovery(R = 50)
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

Lopes, J. E. (2024). *Beta Regression for Interval-Censored
Scale-Derived Outcomes*. MSc Dissertation, PPGMNE/UFPR.

Pinheiro, J. C. and Bates, D. M. (2000). *Mixed-Effects Models in S and
S-PLUS*. Springer. DOI: 10.1007/b98882. Validated online via:
<https://doi.org/10.1007/b98882>.

Rue, H., Martino, S., and Chopin, N. (2009). Approximate Bayesian
inference for latent Gaussian models by using integrated nested Laplace
approximations. *Journal of the Royal Statistical Society: Series B*,
71(2), 319-392. DOI: 10.1111/j.1467-9868.2008.00700.x. Validated online
via: <https://doi.org/10.1111/j.1467-9868.2008.00700.x>.
