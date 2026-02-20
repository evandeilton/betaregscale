# Goodness-of-fit measures

Goodness-of-fit measures

## Usage

``` r
brs_gof(object)
```

## Arguments

- object:

  A fitted `"betaregscale"` object.

## Value

Data frame with logLik, AIC, BIC, and pseudo-R-squared.

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
set.seed(42)
n <- 100
dat <- data.frame(x1 = rnorm(n))
sim <- brs_sim(formula = ~ x1, data = dat, beta = c(0.2, 0.5), phi = 0.3, ncuts = 10)
fit <- brs(y ~ x1, data = sim)
brs_gof(fit)
#>      logLik     AIC      BIC  pseudo_r2
#> 1 -274.5445 555.089 562.9045 0.04425511
# }
```
