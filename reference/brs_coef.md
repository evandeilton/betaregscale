# Internal coefficient table (deprecated, use brs_est() or summary())

Deprecated convenience wrapper. Use
[`brs_est`](https://evandeilton.github.io/betaregscale/reference/brs_est.md)
for coefficient estimates or
[`summary.brs`](https://evandeilton.github.io/betaregscale/reference/summary.brs.md)
for a full model summary.

## Usage

``` r
brs_coef(fit, alpha = 0.05)
```

## Arguments

- fit:

  A fitted `"brs"` object.

- alpha:

  Significance level.

## Value

A list with components `est` (from
[`brs_est`](https://evandeilton.github.io/betaregscale/reference/brs_est.md))
and `gof` (from
[`brs_gof`](https://evandeilton.github.io/betaregscale/reference/brs_gof.md)).

## References

Lopes, J. E. (2023). *Modelos de regressao beta para dados de escala*.
Master's dissertation, Universidade Federal do Parana, Curitiba. URI:
<https://hdl.handle.net/1884/86624>.

Ferrari, S. L. P., and Cribari-Neto, F. (2004). Beta regression for
modelling rates and proportions. *Journal of Applied Statistics*,
**31**(7), 799–815.
[doi:10.1080/0266476042000214501](https://doi.org/10.1080/0266476042000214501)

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

## See also

[`brs_est`](https://evandeilton.github.io/betaregscale/reference/brs_est.md),
[`brs_gof`](https://evandeilton.github.io/betaregscale/reference/brs_gof.md),
[`summary.brs`](https://evandeilton.github.io/betaregscale/reference/summary.brs.md)

## Examples

``` r
# \donttest{
dat <- data.frame(
  y = c(
    0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
    10, 40, 55, 70, 85, 25, 35, 65, 80, 15
  ),
  x1 = rep(c(1, 2), 10)
)
prep <- brs_prep(dat, ncuts = 100)
#> brs_prep: n = 20 | exact = 0, left = 1, right = 1, interval = 18
fit <- brs(y ~ x1, data = prep)
suppressWarnings(brs_coef(fit))  # deprecated; use brs_est()
#> $est
#>      variable   estimate        se    z_value   p_value   ci_lower  ci_upper
#> 1 (Intercept)  0.2551000 0.8643917  0.2951208 0.7679016 -1.4390767 1.9492766
#> 2          x1 -0.2202060 0.5411949 -0.4068885 0.6840899 -1.2809285 0.8405165
#> 3       (phi) -0.3929144 0.2762618 -1.4222538 0.1549526 -0.9343775 0.1485488
#> 
#> $gof
#>      logLik      AIC      BIC   pseudo_r2
#> 1 -92.65206 191.3041 194.2913 0.002947767
#> 
# }
```
