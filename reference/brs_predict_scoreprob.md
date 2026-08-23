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
# \donttest{
dat <- data.frame(
  y = c(
    0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
    10, 40, 55, 70, 85, 25, 35, 65, 80, 15
  ),
  x1 = rep(c(1, 2), 10),
  x2 = rep(c(0, 0, 1, 1), 5)
)
prep <- brs_prep(dat, ncuts = 100)
#> brs_prep: n = 20 | exact = 0, left = 1, right = 1, interval = 18
fit <- brs(y ~ x1, data = prep)
pmat <- brs_predict_scoreprob(fit)
head(pmat[, 1:5])
#>         score_0    score_1    score_2    score_3    score_4
#> [1,] 0.01415137 0.01827159 0.01528005 0.01384045 0.01292155
#> [2,] 0.02412195 0.02639959 0.02075469 0.01816021 0.01654400
#> [3,] 0.01415137 0.01827159 0.01528005 0.01384045 0.01292155
#> [4,] 0.02412195 0.02639959 0.02075469 0.01816021 0.01654400
#> [5,] 0.01415137 0.01827159 0.01528005 0.01384045 0.01292155
#> [6,] 0.02412195 0.02639959 0.02075469 0.01816021 0.01654400
plong <- brs_predict_scoreprob(fit, scores = 0:10, format = "long")
head(plong)
#>   id score       prob
#> 1  1     0 0.01415137
#> 2  1     1 0.01827159
#> 3  1     2 0.01528005
#> 4  1     3 0.01384045
#> 5  1     4 0.01292155
#> 6  1     5 0.01226156
# }
```
