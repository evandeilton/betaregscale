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
#>      variable   estimate        se    z_value      p_value    ci_lower
#> 1 (Intercept) -2.9243597 0.1495297 -19.557049 3.592790e-85 -3.21743259
#> 2          x1  0.1991414 0.1195263   1.666089 9.569577e-02 -0.03512584
#> 3       (phi) -2.8784291 0.2435866 -11.816860 3.194011e-32 -3.35585017
#>     ci_upper
#> 1 -2.6312869
#> 2  0.4334087
#> 3 -2.4010081
# }
```
