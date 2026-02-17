# Marginal effects for brs models

Computes average marginal effects (AME) for numeric covariates in the
mean or precision submodel of a fitted `"brs"` object.

## Usage

``` r
brs_marginaleffects(
  object,
  newdata = NULL,
  model = c("mean", "precision"),
  type = c("response", "link"),
  variables = NULL,
  h = 1e-05,
  interval = TRUE,
  level = 0.95,
  n_sim = 400L,
  seed = NULL,
  keep_draws = FALSE
)
```

## Arguments

- object:

  A fitted `"brs"` object.

- newdata:

  Optional data frame for evaluation; defaults to the data used in
  fitting.

- model:

  Character; `"mean"` (default) or `"precision"`.

- type:

  Character prediction scale: `"response"` (default) or `"link"`.

- variables:

  Optional character vector of covariate names. Defaults to all numeric
  covariates in the selected submodel.

- h:

  Finite-difference step for non-binary numeric covariates.

- interval:

  Logical; compute interval estimates via simulation.

- level:

  Confidence level for interval estimates.

- n_sim:

  Number of parameter draws when `interval = TRUE`.

- seed:

  Optional random seed for reproducibility.

- keep_draws:

  Logical; if `TRUE` and `interval = TRUE`, stores AME simulation draws
  in attribute `"ame_draws"`.

## Value

A data frame with one row per variable and columns: `variable`, `ame`,
`std.error`, `ci.lower`, `ci.upper`, `model`, `type`, and `n`. The
returned object has class `"brs_marginaleffects"` and attributes with
analysis metadata.

## Details

AMEs are computed by finite differences on predictions: \$\$
\mathrm{AME}\_j = \frac{1}{n}\sum\_{i=1}^{n} \frac{\hat{g}\_i(x\_{ij} +
h) - \hat{g}\_i(x\_{ij})}{h}, \$\$ where \\\hat{g}\_i\\ is the selected
prediction scale.

For binary covariates coded as `0/1`, the effect is computed as the
average discrete difference \\\hat{g}(x_j=1)-\hat{g}(x_j=0)\\.

If `interval = TRUE`, uncertainty is approximated by asymptotic
parameter simulation from \\\mathcal{N}(\hat{\theta}, \hat{V})\\.

## References

Hawker, G. A., Mian, S., Kendzerska, T., and French, M. (2011). Measures
of adult pain: Visual Analog Scale for Pain (VAS Pain), Numeric Rating
Scale for Pain (NRS Pain), McGill Pain Questionnaire (MPQ), Short-Form
McGill Pain Questionnaire (SF-MPQ), Chronic Pain Grade Scale (CPGS),
Short Form-36 Bodily Pain Scale (SF-36 BPS), and Measure of Intermittent
and Constant Osteoarthritis Pain (ICOAP). Arthritis Care and Research,
63(S11), S240-S252. doi:10.1002/acr.20543.

Hjermstad, M. J., Fayers, P. M., Haugen, D. F., et al. (2011). Studies
comparing Numerical Rating Scales, Verbal Rating Scales, and Visual
Analogue Scales for assessment of pain intensity in adults: a systematic
literature review. Journal of Pain and Symptom Management, 41(6),
1073-1093. doi:10.1016/j.jpainsymman.2010.08.016.

## Examples

``` r
# \donttest{
set.seed(11)
n <- 150
dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n), z1 = rnorm(n))
sim <- brs_sim(
  formula = ~ x1 + x2 | z1, data = dat,
  beta = c(0.2, -0.5, 0.2), zeta = c(0.2, -0.3),
  ncuts = 100, repar = 2
)
fit <- brs(y ~ x1 + x2 | z1, data = sim, repar = 2)

brs_marginaleffects(fit, model = "mean", type = "response")
#>   variable         ame  std.error     ci.lower    ci.upper model     type   n
#> 1       x1 -0.09214699 0.02300199 -0.138337898 -0.04605006  mean response 150
#> 2       x2  0.04192221 0.02492928 -0.006648427  0.08970841  mean response 150
brs_marginaleffects(fit, model = "precision", type = "link")
#>   variable        ame  std.error   ci.lower   ci.upper     model type   n
#> 1       z1 -0.4099968 0.09462978 -0.5998432 -0.2286651 precision link 150
# }
```
