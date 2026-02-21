# betaregscale <img src="man/figures/logo.png" align="right" height="139" />

[![R-CMD-check](https://github.com/evandeilton/betaregscale/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/evandeilton/betaregscale/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![CRAN status](https://www.r-pkg.org/badges/version/betaregscale)](https://CRAN.R-project.org/package=betaregscale) 
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/betaregscale)](https://cran.r-project.org/package=betaregscale)

## The Methodological Gap

Patient-reported outcome measures (PROMs) on bounded rating scales (e.g., NRS-11, VAS) are widely used in clinical research. Standard analyses treat these bounded, discrete scores as exact continuous values using ordinary least squares (OLS), which ignores natural scale boundaries and misrepresents heteroscedasticity. 

While standard beta regression (e.g., the `betareg` package) respects the $(0,1)$ support, it suffers from two critical limitations:

1. **Interpretability:** It relies on a mean-precision parameterization ($\mu, \phi$) where the precision parameter lacks a direct, clinically intuitive meaning.
2. **Measurement Resolution:** It ignores the discretized nature of rating scales. Selecting "5" on an NRS-11 reflects a measurement interval, not a precise point estimate. Ignoring this leads to underestimated residual variance and biased inference.

## The `betaregscale` Solution

`betaregscale` provides a frequentist, maximum-likelihood framework tailored specifically for bounded scale data. It introduces two major methodological advancements:

1. **Mean-Dispersion (MD) Parameterization:** Reparameterizes the beta distribution strictly in terms of the conditional mean $\mu \in (0,1)$ and a proportional dispersion parameter $\sigma \in (0,1)$.
2. **Interval-Censored Likelihood:** Properly treats discrete scale points as interval-censored data, integrating the beta probability density function over the uncertainty bounds implied by the instrument's resolution.

The package features a compiled **C++ backend** for analytical gradient computation and provides a mixed-effects extension (`brsmm()`) utilizing a multivariate **Laplace approximation** to accommodate repeated measures and clustered data.

---

## Installation

Install the development version directly from GitHub:

```r
# install.packages("remotes")
remotes::install_github("evandeilton/betaregscale")
```

---

## Usage & Technical Workflow

### Data Simulation and Fixed-Effects Modeling

```r
library(betaregscale)

# Simulate interval-censored data with covariate-dependent dispersion
set.seed(42)
d <- data.frame(x1 = rnorm(200), z1 = rnorm(200))
sim <- brs_sim(
  formula = ~ x1 | z1, 
  data = d, 
  beta = c(0.3, 0.4), zeta = c(-2, 0.5), # zeta operates on logit(sigma)
  ncuts = 10, repar = 2
)

# Fit the interval-censored fixed-effects model
fit_fe <- brs(y ~ x1 | z1, data = sim, repar = 2)
summary(fit_fe)
```

### Mixed-Effects Modeling (`brsmm`)

```r
# Simulate clustered data for random intercepts & slopes
d_mm <- data.frame(
  x1 = rnorm(500),
  group = factor(rep(1:50, each = 10))
)
sim_mm <- brs_sim(y ~ x1 | 1, data = d_mm, ncuts = 10, repar = 2)

# Fit mixed-effects model with random intercept
fit_ri <- brsmm(y ~ x1, random = ~ 1 | group, data = sim_mm, repar = 2)

# Fit mixed-effects model with random intercept + slope
fit_rs <- brsmm(y ~ x1, random = ~ 1 + x1 | group, data = sim_mm, repar = 2)

# Likelihood-ratio test for nested model comparison
anova(fit_ri, fit_rs, test = "Chisq")
```

### The Analyst Toolkit & Diagnostics

`betaregscale` is designed for end-to-end clinical reporting, supplying an extensive S3 interface (`print`, `summary`, `coef`, `vcov`, `predict`, `confint`, `ranef`), alongside specialized analyst functions:

```r
# 1. Randomized Quantile Residuals (exact standard normal under correct specification)
res_q <- residuals(fit_ri, type = "rqr")

# 2. Average Marginal Effects (AME) on the response scale
brs_marginaleffects(fit_ri, type = "response")

# 3. Predict probabilities for specific discrete scale categories
brs_predict_scoreprob(fit_ri, scores = 0:10)

# 4. Out-of-sample k-fold cross-validation
brs_cv(y ~ x1, data = sim_mm, k = 5, repeats = 1, repar = 2)

# 5. ggplot2 Diagnostics (Residual vs Fitted, QQ, Scale-Location, Half-normal envelope)
autoplot(fit_ri)
```

---

## Mathematical Framework

### Mean-Dispersion Parameterization
Under the MD parameterization (`repar = 2`), the response $Y_i \sim \text{Beta}(\mu_i, \sigma_i)$ has expected value and variance given by:
$$\text{E}(Y) = \mu, \quad \text{Var}(Y) = \mu(1-\mu)\sigma$$

Both the mean and dispersion can be modeled via link functions ($g$ and $h$) allowing for covariate-dependent heteroscedasticity:
$$g(\mu_i) = x_i^\top \beta, \qquad h(\sigma_i) = z_i^\top \zeta$$

### Interval-Censored Likelihood
Raw scores $y_i^{*} \in \{0, \dots, K\}$ are mapped to the unit interval as $y_i = y_i^{*}/K$, with uncertainty intervals $[l_i, u_i] = [y_i - 1/(2K), y_i + 1/(2K)]$. 


Let $\delta_i \in \{0, 1, 2, 3\}$ indicate the censoring type (exact, left, right, or interval). The complete log-likelihood evaluated in `betaregscale` is:

$$\ell(\theta) = \sum_{i:\delta_i=0} \log f(y_i) + \sum_{i:\delta_i=1} \log F(u_i) + \sum_{i:\delta_i=2} \log\bigl[1 - F(l_i)\bigr] + \sum_{i:\delta_i=3} \log\bigl[F(u_i) - F(l_i)\bigr]$$

where $f(\cdot)$ and $F(\cdot)$ are the beta PDF and CDF.

### Mixed-Effects Extension (`brsmm`)
For grouped or longitudinal data, the mean predictor is extended to include group-specific random effects $\mathbf{b}_j \sim \mathcal{N}(\mathbf{0}, D)$:
$$\eta_{\mu,ij} = x_{ij}^\top \beta + \mathbf{w}_{ij}^\top \mathbf{b}_j$$

The marginal log-likelihood is approximated using a multivariate Laplace approximation:
$$\log L_j(\theta) \approx Q_j(\hat{\mathbf{b}}_j) + \frac{q_b}{2}\log(2\pi) - \frac{1}{2}\log|H_j|$$

---

## Documentation & Learning

To get the most out of `betaregscale`, we highly recommend reading our vignettes. They cover everything from the mathematical foundations to advanced workflows:

* [**Introduction to betaregscale**](https://evandeilton.github.io/betaregscale/articles/brs-intro.html): The basics of interval-censoring and data preparation.
* [**Mixed-Effects Beta Interval Regression**](https://evandeilton.github.io/betaregscale/articles/brs-mm.html): Deep dive into `brsmm()` and random effects.
* [**Analyst Tools**](https://evandeilton.github.io/betaregscale/articles/brs-analyst-tools.html): Using bootstrapping, cross-validation, and marginal effects.
* [**Advanced Workflows**](https://evandeilton.github.io/betaregscale/articles/brs-advanced-workflows.html): A step-by-step guide for high-level users and production environments.

---

## License

MIT &copy; José Evandeilton Lopes
