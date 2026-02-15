# Predict score probabilities from a fitted brs model

Computes predicted probabilities for integer scores on the original
scale \\\\0, 1, \ldots, K\\\\ implied by the fitted beta interval model.

## Usage

``` r
brs_predict_scoreprob(
  object,
  newdata = NULL,
  scores = NULL,
  format = c("matrix", "long"),
  id_col = "id"
)
```

## Arguments

- object:

  A fitted `"brs"` object.

- newdata:

  Optional data frame for prediction.

- scores:

  Optional integer vector of scores to evaluate. Defaults to all scores
  in `0:ncuts`.

- format:

  Output format: `"matrix"` (default) or `"long"`.

- id_col:

  Column name for observation id when `format = "long"`.

## Value

If `format = "matrix"`, a numeric matrix with one row per observation
and one column per requested score.

If `format = "long"`, a long data frame with columns `id_col`, `score`,
and `prob`.

## Details

For a score \\s\\ and \\K =\\ `ncuts`, probabilities are computed as:

- \\P(Y=s)=F(\mathrm{lim}/K)\\ for \\s=0\\,

- \\P(Y=s)=1-F((K-\mathrm{lim})/K)\\ for \\s=K\\,

- \\P(Y=s)=F((s+\mathrm{lim})/K)-F((s-\mathrm{lim})/K)\\ for \\s \in
  \\1,\ldots,K-1\\\\,

where \\F\\ is the beta CDF under the fitted \\(\mu_i,\phi_i)\\.

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
set.seed(33)
dat <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
sim <- brs_sim(
  formula = ~ x1 + x2, data = dat,
  beta = c(0.1, -0.4, 0.3), phi = 0.2, ncuts = 100, repar = 2
)
fit <- brs(y ~ x1 + x2, data = sim, repar = 2)

pmat <- brs_predict_scoreprob(fit)
head(pmat[, 1:5])
#>          score_0     score_1     score_2     score_3     score_4
#> [1,] 0.031430978 0.023248768 0.016145283 0.013227826 0.011526741
#> [2,] 0.074951430 0.040971042 0.026186408 0.020498047 0.017296166
#> [3,] 0.176295577 0.063169280 0.036799106 0.027385744 0.022299473
#> [4,] 0.147774516 0.058628992 0.034844930 0.026215090 0.021510394
#> [5,] 0.004123489 0.004724261 0.003795534 0.003371817 0.003111343
#> [6,] 0.363128382 0.074208336 0.039629538 0.028137509 0.022166138

plong <- brs_predict_scoreprob(fit, scores = 0:10, format = "long")
head(plong)
#>   id score       prob
#> 1  1     0 0.03143098
#> 2  1     1 0.02324877
#> 3  1     2 0.01614528
#> 4  1     3 0.01322783
#> 5  1     4 0.01152674
#> 6  1     5 0.01038186
# }
```
