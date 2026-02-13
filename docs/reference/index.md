# Package index

## Model Fitting

Core functions for fitting beta interval regression models.

- [`betaregscale()`](https://evandeilton.github.io/betaregscale/reference/betaregscale.md)
  : Fit a beta interval regression model
- [`betaregscale_fit()`](https://evandeilton.github.io/betaregscale/reference/betaregscale_fit.md)
  : Fit a fixed-dispersion beta interval regression model
- [`betaregscale_fit_z()`](https://evandeilton.github.io/betaregscale/reference/betaregscale_fit_z.md)
  : Fit a variable-dispersion beta interval regression model

## Log-Likelihood

Log-likelihood evaluation for fixed and variable dispersion.

- [`betaregscale_loglik()`](https://evandeilton.github.io/betaregscale/reference/betaregscale_loglik.md)
  : Log-likelihood for fixed-dispersion beta interval regression
- [`betaregscale_loglik_z()`](https://evandeilton.github.io/betaregscale/reference/betaregscale_loglik_z.md)
  : Log-likelihood for variable-dispersion beta interval regression

## Simulation

Data simulation for Monte Carlo studies.

- [`betaregscale_simulate()`](https://evandeilton.github.io/betaregscale/reference/betaregscale_simulate.md)
  : Simulate data from a fixed-dispersion beta interval model
- [`betaregscale_simulate_z()`](https://evandeilton.github.io/betaregscale/reference/betaregscale_simulate_z.md)
  : Simulate data from a variable-dispersion beta interval model

## S3 Methods

Standard methods for fitted model objects of class betaregscale.

- [`coef(`*`<betaregscale>`*`)`](https://evandeilton.github.io/betaregscale/reference/coef.betaregscale.md)
  : Extract model coefficients
- [`vcov(`*`<betaregscale>`*`)`](https://evandeilton.github.io/betaregscale/reference/vcov.betaregscale.md)
  : Variance-covariance matrix of estimated coefficients
- [`summary(`*`<betaregscale>`*`)`](https://evandeilton.github.io/betaregscale/reference/summary.betaregscale.md)
  : Summarize a fitted model (betareg style)
- [`print(`*`<betaregscale>`*`)`](https://evandeilton.github.io/betaregscale/reference/print.betaregscale.md)
  : Print a fitted model (brief betareg style)
- [`print(`*`<summary.betaregscale>`*`)`](https://evandeilton.github.io/betaregscale/reference/print.summary.betaregscale.md)
  : Print a model summary (betareg style)
- [`logLik(`*`<betaregscale>`*`)`](https://evandeilton.github.io/betaregscale/reference/logLik.betaregscale.md)
  : Extract log-likelihood
- [`AIC(`*`<betaregscale>`*`)`](https://evandeilton.github.io/betaregscale/reference/AIC.betaregscale.md)
  : Akaike information criterion
- [`BIC(`*`<betaregscale>`*`)`](https://evandeilton.github.io/betaregscale/reference/BIC.betaregscale.md)
  : Bayesian information criterion
- [`nobs(`*`<betaregscale>`*`)`](https://evandeilton.github.io/betaregscale/reference/nobs.betaregscale.md)
  : Number of observations
- [`formula(`*`<betaregscale>`*`)`](https://evandeilton.github.io/betaregscale/reference/formula.betaregscale.md)
  : Extract model formula
- [`model.matrix(`*`<betaregscale>`*`)`](https://evandeilton.github.io/betaregscale/reference/model.matrix.betaregscale.md)
  : Extract design matrix
- [`fitted(`*`<betaregscale>`*`)`](https://evandeilton.github.io/betaregscale/reference/fitted.betaregscale.md)
  : Extract fitted values
- [`residuals(`*`<betaregscale>`*`)`](https://evandeilton.github.io/betaregscale/reference/residuals.betaregscale.md)
  : Extract residuals
- [`predict(`*`<betaregscale>`*`)`](https://evandeilton.github.io/betaregscale/reference/predict.betaregscale.md)
  : Predict from a fitted model
- [`confint(`*`<betaregscale>`*`)`](https://evandeilton.github.io/betaregscale/reference/confint.betaregscale.md)
  : Wald confidence intervals
- [`plot(`*`<betaregscale>`*`)`](https://evandeilton.github.io/betaregscale/reference/plot.betaregscale.md)
  : Diagnostic plots for beta interval regression

## Diagnostics and Summaries

Functions for model diagnostics and censoring summaries.

- [`censoring_summary()`](https://evandeilton.github.io/betaregscale/reference/censoring_summary.md)
  : Graphical and tabular censoring summary
- [`est()`](https://evandeilton.github.io/betaregscale/reference/est.md)
  : Coefficient estimates with inference
- [`gof()`](https://evandeilton.github.io/betaregscale/reference/gof.md)
  : Goodness-of-fit measures
- [`hessian_matrix()`](https://evandeilton.github.io/betaregscale/reference/hessian_matrix.md)
  : Extract the Hessian matrix

## Data Preparation

Functions for preparing response variables and reparameterization.

- [`bs_prepare()`](https://evandeilton.github.io/betaregscale/reference/bs_prepare.md)
  : Pre-process analyst data for beta interval regression
- [`check_response()`](https://evandeilton.github.io/betaregscale/reference/check_response.md)
  : Transform and validate a scale-derived response variable
- [`beta_reparam()`](https://evandeilton.github.io/betaregscale/reference/beta_reparam.md)
  : Reparameterize (mu, phi) into beta shape parameters
