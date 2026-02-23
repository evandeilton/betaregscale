# betaregscale <img src="man/figures/logo.png" align="right" height="139" />

[![R-CMD-check](https://github.com/evandeilton/betaregscale/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/evandeilton/betaregscale/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![CRAN status](https://www.r-pkg.org/badges/version/betaregscale)](https://CRAN.R-project.org/package=betaregscale)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/betaregscale)](https://cran.r-project.org/package=betaregscale)

## The Methodological Gap

Patient-reported outcome measures (PROMs) on bounded rating scales — such as the
**Numerical Rating Scale (NRS-11, NRS-21, NRS-101)**, Visual Analogue Scale (VAS),
or Likert-type instruments — are ubiquitous in clinical research. Standard analyses
treat these bounded, discrete scores as exact continuous values via ordinary least
squares (OLS), silently violating two core assumptions: natural scale boundaries
are ignored, and the inherent heteroscedasticity of bounded responses is
misrepresented.

While standard beta regression (e.g., the `betareg` package) respects the $(0,1)$
support, it suffers from two critical limitations:

1.  **Interpretability:** It relies on a mean-precision parameterization $(\mu,\phi)$
    where the precision parameter $\phi$ lacks a direct, clinically intuitive meaning.
2.  **Measurement Resolution:** It ignores the discretized nature of rating scales.
    Selecting "52" on an NRS-101 reflects a measurement interval $[0.515, 0.525]$,
    not a precise continuous value. Ignoring this leads to underestimated residual
    variance and biased inference.

## The `betaregscale` Solution

`betaregscale` provides a frequentist, maximum-likelihood framework tailored
specifically for bounded scale data. It introduces two major methodological
advancements:

1.  **Mean-Dispersion (MD) Parameterization:** Reparameterizes the beta distribution
    in terms of the conditional mean $\mu \in (0,1)$ and a proportional dispersion
    parameter $\sigma \in (0,1)$ — both directly interpretable and modelable via
    covariates.

2.  **Interval-Censored Likelihood:** Properly treats each discrete scale point as
    interval-censored data, integrating the beta PDF over the uncertainty bounds
    implied by the instrument's resolution. A score of $y^*$ on a $K$-point scale
    is treated as $[y^*/K - 1/(2K),\; y^*/K + 1/(2K)]$.

The package features a compiled **C++ backend** for analytical gradient computation,
and provides a mixed-effects extension (`brsmm()`) via multivariate **Laplace
approximation** for repeated measures and multi-centre data.

---

## Installation

```r
# Stable version from CRAN
install.packages("betaregscale")

# Development version from GitHub
# install.packages("remotes")
remotes::install_github("evandeilton/betaregscale")
```

---

## Usage Examples

```r
# =============================================================================
# betaregscale — Clinical Workflow: NRS-101 Pain Score Modelling
# =============================================================================
# Scenario: multi-centre analgesic RCT (500 patients, 20 clinics).
# Outcome:  post-treatment pain rated on NRS-101 (0 = no pain, 100 = worst).
# Goal:     illustrate the full betaregscale pipeline — simulation, fixed-
#           effects modelling, mixed-effects modelling, and post-estimation
#           diagnostics.
# =============================================================================

library(betaregscale)
library(ggplot2)


# -----------------------------------------------------------------------------
# 1. CLINICAL PREDICTORS
# -----------------------------------------------------------------------------
# Simulate a realistic patient-level covariate matrix:
#   trt       — treatment arm (0 = placebo, 1 = active compound; 50/50)
#   age_std   — age in years (25–75), standardised (mean 0, sd 1)
#   sex_f     — sex (1 = female, 0 = male; ~55 % female)
#   bpain_std — baseline NRS-101 pain (30–90), standardised
#   clinic    — clinic ID (20 centres, 25 patients each)

set.seed(2026)
n_patients <- 1000

patient_data <- data.frame(
  trt       = rbinom(n_patients, 1, prob = 0.50),
  age_std   = scale(runif(n_patients, 25, 75))[, 1],
  sex_f     = rbinom(n_patients, 1, prob = 0.55),
  bpain_std = scale(runif(n_patients, 30, 90))[, 1],
  clinic    = factor(rep(1:20, each = 25))
)


# -----------------------------------------------------------------------------
# 2. DATA SIMULATION WITH brs_sim()
# -----------------------------------------------------------------------------
# Generate NRS-101 scores from the true beta-regression model.
# brs_sim() returns a data frame already containing the interval-censored
# columns (left, right, yt, delta) required by brs().
#
# True parameters — mean sub-model (logit link on mu):
#   Intercept =  0.20  → average pain ~55/100 in the placebo arm
#   trt       = -0.55  → active drug reduces pain by ~13 NRS points
#   age_std   =  0.15  → older patients report slightly higher pain
#   sex_f     =  0.10  → females report slightly higher pain
#   bpain_std =  0.60  → strong tracking from baseline (dominant predictor)
#
# True parameters — dispersion sub-model (logit link on sigma):
#   Intercept = -1.60  → moderate dispersion (~17 %)
#   trt       = -0.40  → active arm shows more homogeneous responses
#   bpain_std =  0.20  → high baseline pain leads to wider response spread

sim_nrs101 <- brs_sim(
  formula = ~ trt + age_std + sex_f + bpain_std | trt + bpain_std,
  data    = patient_data,
  beta    = c(0.20, -0.55, 0.15, 0.10, 0.60),
  zeta    = c(-1.60, -0.40, 0.20),
  ncuts   = 100,   # NRS-101: integer scores 0–100 → ncuts = 100
  repar   = 2      # Mean-Dispersion (MD) parameterisation
)

# Attach clinic for the mixed-effects section
sim_nrs101$clinic <- patient_data$clinic


# -----------------------------------------------------------------------------
# 3. FIXED-EFFECTS MODELS
# -----------------------------------------------------------------------------

# 3a. Full model — covariate-dependent dispersion (correctly specified)
#     Mean:       mu    ~ trt + age_std + sex_f + bpain_std   (logit link)
#     Dispersion: sigma ~ trt + bpain_std                     (logit link)
fit_full <- brs(
  formula  = y ~ trt + age_std + sex_f + bpain_std | trt + bpain_std,
  data     = sim_nrs101,
  repar    = 2,
  link     = "logit",
  link_phi = "logit"
)

summary(fit_full)

# 3b. Parsimonious model — constant (fixed) dispersion
#     Used as baseline for model comparison.
fit_fixed <- brs(
  formula = y ~ trt + age_std + sex_f + bpain_std,
  data    = sim_nrs101,
  repar   = 2
)

# 3c. Model comparison table (AIC / BIC / log-likelihood / pseudo-R²)
#     Variable dispersion should win if the DGP is correctly recovered.
brs_table(
  "Fixed dispersion"    = fit_fixed,
  "Variable dispersion" = fit_full,
  sort_by = "AIC"
)


# -----------------------------------------------------------------------------
# 4. MIXED-EFFECTS MODELS WITH brsmm()
# -----------------------------------------------------------------------------
# Clinics introduce between-centre heterogeneity; ignoring this structure
# inflates Type I error. brsmm() adds Gaussian random effects via Laplace
# approximation.

# 4a. Random-intercept model — clinics differ in average pain level
fit_ri <- brsmm(
  formula = y ~ trt + age_std + sex_f + bpain_std | trt + bpain_std,
  random  = ~ 1 | clinic,
  data    = sim_nrs101,
  repar   = 2
)

# 4b. Random intercept + random slope — treatment efficacy may vary by clinic
fit_rs <- brsmm(
  formula = y ~ trt + age_std + sex_f + bpain_std | trt + bpain_std,
  random  = ~ 1 + trt | clinic,
  data    = sim_nrs101,
  repar   = 2
)

# 4c. Likelihood-ratio test: is a random treatment slope warranted?
#     p < 0.05 → clinics differ significantly in treatment response.
anova(fit_ri, fit_rs, test = "Chisq")

# 4d. Clinic-level BLUPs (posterior modes of random effects)
ranef(fit_ri)


# -----------------------------------------------------------------------------
# 5. POST-ESTIMATION: RESIDUALS & DIAGNOSTICS
# -----------------------------------------------------------------------------

# 5a. Randomized quantile residuals (RQRs)
#     Under the correctly specified model, RQRs ~ N(0,1) exactly.
#     Departure from the diagonal in the Q-Q plot signals misspecification.
rqr <- residuals(fit_full, type = "rqr")
qqnorm(rqr); qqline(rqr, col = "steelblue", lwd = 2)

# 5b. Four-panel ggplot2 diagnostic grid:
#     residuals vs fitted | Q-Q | scale-location | half-normal envelope
autoplot(fit_full)

# 5c. Calibration plot: predicted mu (decile bins) vs. empirical mean
autoplot(fit_full, type = "calibration")

# 5d. Predicted score distribution over the full NRS-101 grid (0–100)
autoplot(fit_full, type = "score_dist", scores = 0:100)


# -----------------------------------------------------------------------------
# 6. AVERAGE MARGINAL EFFECTS (AME)
# -----------------------------------------------------------------------------
# AMEs quantify the average unit change in the predicted mean (response scale,
# mu ∈ [0,1]) for a one-unit change in each predictor.
# Multiply by 100 to interpret directly in NRS-101 points.

set.seed(2026)
ame <- brs_marginaleffects(
  fit_full,
  model    = "mean",
  type     = "response",
  interval = TRUE,
  n_sim    = 200
)
ame


# -----------------------------------------------------------------------------
# 7. SCORE-LEVEL PREDICTED PROBABILITIES
# -----------------------------------------------------------------------------
# For each patient, compute the probability mass assigned to each NRS-101 bin.
# Here we evaluate a coarse 11-point grid (multiples of 10) for readability;
# pass scores = 0:100 for the full resolution.

prob_scores <- brs_predict_scoreprob(
  fit_full,
  scores = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
)

# 500 patients × 11 score bins
dim(prob_scores)

# Probability profile for the first 10 patients
round(head(prob_scores, 10), 3)


# -----------------------------------------------------------------------------
# 8. CROSS-VALIDATION
# -----------------------------------------------------------------------------
# Repeated 5-fold CV evaluates out-of-sample predictive performance.
# Reports: log-score, RMSE_yt, and MAE_yt — averaged across folds and repeats.

set.seed(2203)
cv_res <- brs_cv(
  formula = y ~ trt + age_std + sex_f + bpain_std | trt + bpain_std,
  data    = sim_nrs101,
  k       = 5,
  repeats = 3,
  repar   = 2
)
cv_res


# -----------------------------------------------------------------------------
# 9. CONFIDENCE INTERVALS & BACK-TRANSFORMATION TO NRS-101
# -----------------------------------------------------------------------------

# 9a. 95 % Wald CIs for all fixed-effect parameters
confint(fit_full)

# 9b. Patient-level predicted means on the mu ∈ [0,1] scale
mu_hat <- predict(fit_full, type = "response")

# 9c. Back-transform to NRS-101 units (multiply by 100)
nrs101_hat <- mu_hat * 100
summary(nrs101_hat)
```

---

## Mathematical Framework

### Mean-Dispersion Parameterization

Under `repar = 2`, the response $Y_i \sim \text{Beta}(\mu_i, \sigma_i)$ has:

$$
\text{E}(Y_i) = \mu_i, \qquad \text{Var}(Y_i) = \mu_i(1-\mu_i)\sigma_i
$$

Both submodels are linked to linear predictors via monotone link functions
$g$ and $h$:

$$
g(\mu_i) = \mathbf{x}_i^\top \boldsymbol{\beta}, \qquad h(\sigma_i) = \mathbf{z}_i^\top \boldsymbol{\zeta}
$$

The default link for both is logit, ensuring $\mu_i, \sigma_i \in (0,1)$.

### Interval-Censored Likelihood

A raw NRS-101 score $y_i^* \in \{0, 1, \ldots, 100\}$ is mapped to
$y_i = y_i^*/100$ with interval $[l_i, u_i] = [y_i - 0.005,\, y_i + 0.005]$.

Let $\delta_i \in \{0,1,2,3\}$ encode the censoring type. The complete
log-likelihood is:

$$
\ell(\boldsymbol{\theta}) =
  \sum_{i:\,\delta_i=0} \log f(y_i)
+ \sum_{i:\,\delta_i=1} \log F(u_i)
+ \sum_{i:\,\delta_i=2} \log\bigl[1 - F(l_i)\bigr]
+ \sum_{i:\,\delta_i=3} \log\bigl[F(u_i) - F(l_i)\bigr]
$$

where $f(\cdot)$ and $F(\cdot)$ are the Beta PDF and CDF evaluated at the
shape parameters derived from $(\mu_i, \sigma_i)$.

### Mixed-Effects Extension (`brsmm`)

For patient $i$ in clinic $j$, the mean predictor is extended by
clinic-specific random effects $\mathbf{b}_j \sim \mathcal{N}(\mathbf{0}, D)$:

$$
\eta_{\mu,ij} = \mathbf{x}_{ij}^\top \boldsymbol{\beta} + \mathbf{w}_{ij}^\top \mathbf{b}_j
$$

The intractable group marginal likelihood is approximated via Laplace:

$$
\log L_j(\boldsymbol{\theta}) \approx
  Q_j(\hat{\mathbf{b}}_j)
+ \frac{q_b}{2}\log(2\pi)
- \frac{1}{2}\log|H_j|
$$

where $q_b$ is the random-effects dimension, $\hat{\mathbf{b}}_j$ is the
posterior mode, and $H_j = -\nabla^2 Q_j(\hat{\mathbf{b}}_j)$.

---

## S3 Interface Summary

| Generic | `brs` | `brsmm` |
|---|:---:|:---:|
| `print()` | ✓ | ✓ |
| `summary()` | ✓ | ✓ |
| `coef()` | ✓ | ✓ |
| `vcov()` | ✓ | ✓ |
| `confint()` | ✓ | ✓ |
| `predict()` | ✓ | ✓ |
| `residuals()` (RQR, response, pearson) | ✓ | ✓ |
| `ranef()` | — | ✓ |
| `anova()` | ✓ | ✓ |
| `autoplot()` | ✓ | — |
| `brs_marginaleffects()` | ✓ | — |
| `brs_predict_scoreprob()` | ✓ | — |
| `brs_cv()` | ✓ | — |
| `brs_table()` | ✓ | ✓ |

---

## Documentation & Vignettes

| Vignette | Content |
|---|---|
| [**Introduction to betaregscale**](https://evandeilton.github.io/betaregscale/articles/brs-intro.html) | Censoring types, `brs_prep()` modes, first fixed-effects model |
| [**Mixed-Effects Beta Interval Regression**](https://evandeilton.github.io/betaregscale/articles/brs-mm.html) | `brsmm()` mathematics, parameter recovery, BLUP extraction |
| [**Analyst Tools**](https://evandeilton.github.io/betaregscale/articles/brs-analyst-tools.html) | AME, score probabilities, cross-validation, `autoplot()` |
| [**Advanced Workflows**](https://evandeilton.github.io/betaregscale/articles/brs-advanced-workflows.html) | Production pipelines, model selection, sensitivity analysis |

---

## Citation

If you use `betaregscale` in a publication, please cite:

```r
citation("betaregscale")
```

---

## License

MIT © José Evandeilton Lopes & Wagner Hugo Bonat
