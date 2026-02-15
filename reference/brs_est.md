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
#>      variable   estimate        se    z_value      p_value     ci_lower
#> 1 (Intercept) -3.0886607 0.1787335 -17.280815 6.561480e-67 -3.438971974
#> 2          x1  0.3326557 0.1729164   1.923794 5.438036e-02 -0.006254317
#> 3       (phi) -2.7342024 0.2646231 -10.332439 5.026571e-25 -3.252854133
#>     ci_upper
#> 1 -2.7383495
#> 2  0.6715657
#> 3 -2.2155506
# }
```
