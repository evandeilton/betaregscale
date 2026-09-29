# betaregscale: Beta Regression for Interval-Censored Scale-Derived Outcomes

Maximum-likelihood beta regression for scores recorded on a bounded
scale \\\\0, 1, \ldots, K\\\\ (pain rating scales, Likert-type items,
ratings), where \\K =\\ `ncuts` is the maximum score; a scale that
starts at 1 is shifted to start at 0. A score is read as a coarsened
observation of a latent \\Y \in (0, 1)\\ with a beta distribution: each
score maps to a cell \\\[l_i, u_i\]\\ of \\(0, 1)\\ and contributes the
beta probability of that cell to the likelihood. This is the interval
beta regression model of Lopes (2023). The package fits fixed- and
variable-dispersion models
([`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md))
and mixed models with random intercepts and slopes
([`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md)),
and provides simulation, score probabilities, bootstrap,
cross-validation and marginal effects.

## Model

For observation \\i\\, \$\$Y_i \sim \mathrm{Beta}(a_i, b_i), \qquad
g_1(\mu_i) = x_i^\top \beta, \qquad g_2(\phi_i) = z_i^\top \gamma,\$\$
where \\(a_i, b_i)\\ follow from \\(\mu_i, \phi_i)\\ under one of three
parameterisations
([`brs_repar`](https://evandeilton.github.io/betaregscale/reference/brs_repar.md)).
The default, `repar = 2`, uses the mean \\\mu\\ and the dispersion
\\\phi = 1/(1 + a + b)\\, written \\\sigma\\ in Lopes (2023). The score
\\s_i\\ gives the cell \\\[l_i, u_i\]\\ and the censoring type
\\\delta_i\\
([`brs_check`](https://evandeilton.github.io/betaregscale/reference/brs_check.md)),
and the log-likelihood is \$\$\ell(\beta, \gamma) = \sum\_{\delta_i = 0}
\log f(y_i) + \sum\_{\delta_i = 1} \log F(u_i) + \sum\_{\delta_i = 2}
\log\\1 - F(l_i)\\ + \sum\_{\delta_i = 3} \log\\F(u_i) - F(l_i)\\,\$\$
with \\f\\ and \\F\\ the beta density and distribution function. Scores
0 and \\K\\ are left- and right-censored (\\\delta = 1, 2\\), other
scores interval-censored (\\\delta = 3\\); values already in \\(0, 1)\\
are exact (\\\delta = 0\\).

## Main functions

- [`brs_check`](https://evandeilton.github.io/betaregscale/reference/brs_check.md),
  [`brs_prep`](https://evandeilton.github.io/betaregscale/reference/brs_prep.md):

  Scores (or analyst-supplied intervals) to cells and censoring types.

- [`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md):

  Model fit;
  [`brs_fit_fixed`](https://evandeilton.github.io/betaregscale/reference/brs_fit_fixed.md)
  and
  [`brs_fit_var`](https://evandeilton.github.io/betaregscale/reference/brs_fit_var.md)
  are the fixed- and variable-dispersion workers.

- [`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md):

  Mixed model with Gaussian random intercepts and slopes (Laplace,
  adaptive Gauss-Hermite or quasi-Monte Carlo integration).

- [`brs_sim`](https://evandeilton.github.io/betaregscale/reference/brs_sim.md):

  Simulate scores from the model.

- [`brs_predict_scoreprob`](https://evandeilton.github.io/betaregscale/reference/brs_predict_scoreprob.md):

  Predicted probability of each score.

- [`brs_bootstrap`](https://evandeilton.github.io/betaregscale/reference/brs_bootstrap.md),
  [`brs_cv`](https://evandeilton.github.io/betaregscale/reference/brs_cv.md),
  [`brs_marginaleffects`](https://evandeilton.github.io/betaregscale/reference/brs_marginaleffects.md),
  [`brs_table`](https://evandeilton.github.io/betaregscale/reference/brs_table.md):

  Parametric bootstrap, cross-validation, average marginal effects and
  model comparison tables.

Fits of class `"brs"` and `"brsmm"` have methods for `print`, `summary`,
`coef`, `vcov`, `confint`, `logLik`, `AIC`, `BIC`, `nobs`, `anova`,
`fitted`, `residuals`, `predict`, `plot` and `autoplot`. As in betareg,
[`coef()`](https://rdrr.io/r/stats/coef.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) take
`model = c("full", "mean", "precision")`. Every fit checks its gradient,
Hessian and likelihood clamps and warns in one line when something is
wrong ('Fit diagnostics' in
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md)).

## Relation to other approaches

The model is the interval beta regression (model M2) of Lopes (2023).
That dissertation compares it with beta regression on the rescaled
scores (betareg, M1) and with the quasi-beta regression of Bonat et al.
(2019), fitted with mcglm (Bonat and Jørgensen, 2016), which specifies
only the mean and the variance (M3). In that study's simulations M3 had
Wald coverage closest to 95%, M1 and M2 fell below 90% in several
scenarios with more than 500 observations and dispersion above 0.2, and
M3 tended to underestimate the covariate effects; on its knee-surgery
data M1 and M2 gave similar time effects. These are results of that
study and of the implementation used then, not properties guaranteed by
the package (a Monte Carlo study of the current code is in
[`vignette("brs-advanced-workflows")`](https://evandeilton.github.io/betaregscale/articles/brs-advanced-workflows.md)).
The interval model is useful when the coarsening of the scale is part of
the question: it gives score probabilities, treats the borders of the
scale as censoring and supports likelihood-ratio tests.

## References

Lopes, J. E. (2023). *Modelos de regressao beta para dados de escala*.
Master's dissertation, Universidade Federal do Parana, Curitiba. URI:
https://hdl.handle.net/1884/86624.

Ferrari, S. L. P., and Cribari-Neto, F. (2004). Beta regression for
modelling rates and proportions. *Journal of Applied Statistics*,
**31**(7), 799–815.
[doi:10.1080/0266476042000214501](https://doi.org/10.1080/0266476042000214501)

Bonat, W. H., Petterle, R. R., Hinde, J., and Demétrio, C. G. B. (2019).
Flexible quasi-beta regression models for continuous bounded data.
*Statistical Modelling*, **19**(6), 617–633.

Bonat, W. H., and Jørgensen, B. (2016). Multivariate covariance
generalized linear models. *Journal of the Royal Statistical Society:
Series C (Applied Statistics)*, **65**(5), 649–675.

## See also

[`vignette("brs-intro", package = "betaregscale")`](https://evandeilton.github.io/betaregscale/articles/brs-intro.md)

## Author

**Maintainer**: José Evandeilton Lopes <evandeilton@gmail.com>
([ORCID](https://orcid.org/0009-0007-5887-4084))

Authors:

- José Evandeilton Lopes <evandeilton@gmail.com>
  ([ORCID](https://orcid.org/0009-0007-5887-4084))

- Wagner Hugo Bonat ([ORCID](https://orcid.org/0000-0002-0349-7054))
