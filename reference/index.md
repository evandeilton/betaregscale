# Package index

## Model Fitting

Core functions for fitting beta interval regression models.

- [`brs()`](https://evandeilton.github.io/betaregscale/reference/brs.md)
  : Fit a beta interval regression model
- [`brs_fit_fixed()`](https://evandeilton.github.io/betaregscale/reference/brs_fit_fixed.md)
  : Fit a fixed-dispersion beta interval regression model
- [`brs_fit_var()`](https://evandeilton.github.io/betaregscale/reference/brs_fit_var.md)
  : Fit a variable-dispersion beta interval regression model
- [`brs_coef()`](https://evandeilton.github.io/betaregscale/reference/brs_coef.md)
  : Internal coefficient table (deprecated, use brs_est() or summary())

## Simulation

Data simulation for Monte Carlo studies.

- [`brs_sim()`](https://evandeilton.github.io/betaregscale/reference/brs_sim.md)
  : Simulate data from beta interval models

## S3 Methods

Standard methods for fitted model objects of class brs.

- [`coef(`*`<brs>`*`)`](https://evandeilton.github.io/betaregscale/reference/coef.brs.md)
  : Extract model coefficients
- [`vcov(`*`<brs>`*`)`](https://evandeilton.github.io/betaregscale/reference/vcov.brs.md)
  : Variance-covariance matrix of estimated coefficients
- [`summary(`*`<brs>`*`)`](https://evandeilton.github.io/betaregscale/reference/summary.brs.md)
  : Summarize a fitted model (betareg style)
- [`print(`*`<brs>`*`)`](https://evandeilton.github.io/betaregscale/reference/print.brs.md)
  : Print a fitted model (brief betareg style)
- [`print(`*`<summary.brs>`*`)`](https://evandeilton.github.io/betaregscale/reference/print.summary.brs.md)
  : Print a model summary (betareg style)
- [`logLik(`*`<brs>`*`)`](https://evandeilton.github.io/betaregscale/reference/logLik.brs.md)
  : Extract log-likelihood
- [`AIC(`*`<brs>`*`)`](https://evandeilton.github.io/betaregscale/reference/AIC.brs.md)
  : Akaike information criterion
- [`BIC(`*`<brs>`*`)`](https://evandeilton.github.io/betaregscale/reference/BIC.brs.md)
  : Bayesian information criterion
- [`nobs(`*`<brs>`*`)`](https://evandeilton.github.io/betaregscale/reference/nobs.brs.md)
  : Number of observations
- [`formula(`*`<brs>`*`)`](https://evandeilton.github.io/betaregscale/reference/formula.brs.md)
  : Extract model formula
- [`model.matrix(`*`<brs>`*`)`](https://evandeilton.github.io/betaregscale/reference/model.matrix.brs.md)
  : Extract design matrix
- [`fitted(`*`<brs>`*`)`](https://evandeilton.github.io/betaregscale/reference/fitted.brs.md)
  : Extract fitted values
- [`residuals(`*`<brs>`*`)`](https://evandeilton.github.io/betaregscale/reference/residuals.brs.md)
  : Extract residuals
- [`predict(`*`<brs>`*`)`](https://evandeilton.github.io/betaregscale/reference/predict.brs.md)
  : Predict from a fitted model
- [`confint(`*`<brs>`*`)`](https://evandeilton.github.io/betaregscale/reference/confint.brs.md)
  : Wald confidence intervals
- [`plot(`*`<brs>`*`)`](https://evandeilton.github.io/betaregscale/reference/plot.brs.md)
  : Diagnostic plots for beta interval regression

## Diagnostics and Summaries

Functions for model diagnostics and censoring summaries.

- [`brs_cens()`](https://evandeilton.github.io/betaregscale/reference/brs_cens.md)
  : Graphical and tabular censoring summary
- [`brs_est()`](https://evandeilton.github.io/betaregscale/reference/brs_est.md)
  : Coefficient estimates with inference
- [`brs_gof()`](https://evandeilton.github.io/betaregscale/reference/brs_gof.md)
  : Goodness-of-fit measures
- [`brs_hessian()`](https://evandeilton.github.io/betaregscale/reference/brs_hessian.md)
  : Extract the Hessian matrix
- [`autoplot.brs()`](https://evandeilton.github.io/betaregscale/reference/autoplot.brs.md)
  : ggplot2 autoplot for brs models

## Analyst Tools

Post-estimation tables, effects, score probabilities, and validation.

- [`brs_table()`](https://evandeilton.github.io/betaregscale/reference/brs_table.md)
  : Compare fitted brs models in a single table
- [`brs_marginaleffects()`](https://evandeilton.github.io/betaregscale/reference/brs_marginaleffects.md)
  : Marginal effects for brs models
- [`brs_predict_scoreprob()`](https://evandeilton.github.io/betaregscale/reference/brs_predict_scoreprob.md)
  : Predict score probabilities from a fitted brs model
- [`brs_cv()`](https://evandeilton.github.io/betaregscale/reference/brs_cv.md)
  : K-fold cross-validation for brs models

## Data Preparation

Functions for preparing response data and beta reparameterization.

- [`brs_prep()`](https://evandeilton.github.io/betaregscale/reference/brs_prep.md)
  : Pre-process analyst data for beta interval regression
- [`brs_check()`](https://evandeilton.github.io/betaregscale/reference/brs_check.md)
  : Transform and validate a scale-derived response variable
- [`brs_repar()`](https://evandeilton.github.io/betaregscale/reference/brs_repar.md)
  : Reparameterize (mu, phi) into beta shape parameters
