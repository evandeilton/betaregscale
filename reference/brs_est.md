# Coefficient estimates with inference

Coefficient estimates with inference

## Usage

``` r
brs_est(object, alpha = 0.05)
```

## Arguments

- object:

  A fitted `"betaregscale"` object.

- alpha:

  Significance level (default 0.05).

## Value

Data frame of estimates, standard errors, z-values, and p-values.

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

[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md),
[`brs_gof`](https://evandeilton.github.io/betaregscale/reference/brs_gof.md),
[`brs_hessian`](https://evandeilton.github.io/betaregscale/reference/brs_hessian.md),
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
brs_est(fit)
#>      variable   estimate        se    z_value   p_value   ci_lower  ci_upper
#> 1 (Intercept)  0.2550945 0.8643907  0.2951148 0.7679062 -1.4390802 1.9492691
#> 2          x1 -0.2202075 0.5411943 -0.4068917 0.6840875 -1.2809288 0.8405139
#> 3       (phi) -0.3929198 0.2762612 -1.4222765 0.1549460 -0.9343818 0.1485422
# }
```
