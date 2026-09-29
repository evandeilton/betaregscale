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

- keep_draws:

  Logical; if `TRUE` and `interval = TRUE`, stores AME simulation draws
  in attribute `"ame_draws"`.

## Value

A data frame with one row per variable and columns: `variable`, `ame`,
`std.error`, `ci.lower`, `ci.upper`, `model`, `type`, and `n`. The
returned object has class `"brs_marginaleffects"` and attributes with
analysis metadata.

## Details

AMEs for a numeric covariate are computed by a central difference on
predictions, with the step size scaled by the covariate's standard
deviation: \$\$ \mathrm{AME}\_j = \frac{1}{n}\sum\_{i=1}^{n}
\frac{\hat{g}\_i(x\_{ij} + h_j) - \hat{g}\_i(x\_{ij} - h_j)}{2 h_j},
\qquad h_j = h \cdot \max(\mathrm{sd}(x\_{\cdot j}), 1), \$\$ where
\\\hat{g}\_i\\ is the selected prediction scale and `h` is the (small)
base step supplied by the caller.

For binary covariates coded as `0/1`, the effect is computed as the
average discrete difference \\\hat{g}(x_j=1)-\hat{g}(x_j=0)\\.

If `interval = TRUE`, uncertainty is approximated by asymptotic
parameter simulation from \\\mathcal{N}(\hat{\theta}, \hat{V})\\.

## References

Lopes, J. E. (2023). *Modelos de regressao beta para dados de escala*.
Master's dissertation, Universidade Federal do Parana, Curitiba. URI:
https://hdl.handle.net/1884/86624.

Hawker, G. A., Mian, S., Kendzerska, T., and French, M. (2011). Measures
of adult pain: Visual Analog Scale for Pain (VAS Pain), Numeric Rating
Scale for Pain (NRS Pain), McGill Pain Questionnaire (MPQ), Short-Form
McGill Pain Questionnaire (SF-MPQ), Chronic Pain Grade Scale (CPGS),
Short Form-36 Bodily Pain Scale (SF-36 BPS), and Measure of Intermittent
and Constant Osteoarthritis Pain (ICOAP). Arthritis Care and Research,
63(S11), S240-S252.
[doi:10.1002/acr.20543](https://doi.org/10.1002/acr.20543)

Hjermstad, M. J., Fayers, P. M., Haugen, D. F., et al. (2011). Studies
comparing Numerical Rating Scales, Verbal Rating Scales, and Visual
Analogue Scales for assessment of pain intensity in adults: a systematic
literature review. Journal of Pain and Symptom Management, 41(6),
1073-1093.
[doi:10.1016/j.jpainsymman.2010.08.016](https://doi.org/10.1016/j.jpainsymman.2010.08.016)

## Examples

``` r
# Synthetic NRS-11 scores with a 0/1 time dummy and a numeric covariate (age)
set.seed(3)
nrs <- data.frame(t12 = rep(0:1, each = 60), age = round(rnorm(120, 40, 12)))
shp <- brs_repar(mu = plogis(-1.3 + 0.75 * nrs$t12 + 0.03 * (nrs$age - 40)),
                 phi = 0.3, repar = 2)
nrs$y <- round(10 * rbeta(nrow(nrs), shp$shape1, shp$shape2))
fit <- brs(y ~ t12 + age, data = nrs, ncuts = 10)

# Average marginal effects on E[Y] (0/1 dummy: discrete change; age: slope),
# with intervals from 50 draws of the estimates (use n_sim >= 400 in practice)
set.seed(6)
brs_marginaleffects(fit, n_sim = 50)
#>   variable         ame   std.error     ci.lower    ci.upper model     type   n
#> 1      t12 0.224844538 0.043324091  0.148876732 0.303261197  mean response 120
#> 2      age 0.002873628 0.002248534 -0.001835917 0.006042961  mean response 120

# On the 0-10 score scale (mid cells): multiply by K = 10
set.seed(6)
ame <- brs_marginaleffects(fit, n_sim = 50)
transform(ame[, c("variable", "ame", "ci.lower", "ci.upper")],
          ame = 10 * ame, ci.lower = 10 * ci.lower, ci.upper = 10 * ci.upper)
#>   variable        ame    ci.lower   ci.upper
#> 1      t12 2.24844538  1.48876732 3.03261197
#> 2      age 0.02873628 -0.01835917 0.06042961
```
