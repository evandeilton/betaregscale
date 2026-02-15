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
sim <- brs_sim(
  formula = ~x1, data = data.frame(x1 = rnorm(50)),
  beta = c(0, 0.5), phi = 0.1, ncuts = 10, repar = 2
)
fit <- brs(y ~ x1, data = sim, repar = 2)
brs_est(fit)
#>      variable   estimate        se   z_value      p_value    ci_lower
#> 1 (Intercept) -2.8611765 0.1361136 -21.02050 4.258938e-98 -3.12795428
#> 2          x1  0.2317141 0.1088125   2.12948 3.321453e-02  0.01844551
#> 3       (phi) -3.0117569 0.2372581 -12.69401 6.383324e-37 -3.47677434
#>     ci_upper
#> 1 -2.5943987
#> 2  0.4449827
#> 3 -2.5467395
# }
```
