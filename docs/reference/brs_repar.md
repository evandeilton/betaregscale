# Reparameterize (mu, phi) into beta shape parameters

Converts a mean–dispersion pair \\(\mu, \phi)\\ to the shape parameters
\\(a, b)\\ of the beta distribution under one of three
reparameterization schemes.

## Usage

``` r
brs_repar(mu, phi, repar = 2L)
```

## Arguments

- mu:

  Numeric vector of mean values in \\(0, 1)\\.

- phi:

  Numeric vector (or scalar) of dispersion values.

- repar:

  Integer (0, 1, or 2) selecting the scheme.

## Value

A `data.frame` with columns `shape1` and `shape2`.

## Details

The three schemes are:

- `repar = 0`:

  Direct: \\a = \mu,\\ b = \phi\\.

- `repar = 1`:

  Ferrari–Cribari-Neto: \\a = \mu\phi,\\ b = (1 - \mu)\phi\\, where
  \\\phi\\ acts as a precision parameter.

- `repar = 2`:

  Mean–variance: \\a = \mu(1-\phi)/\phi,\\ b = (1-\mu)(1-\phi)/\phi\\,
  where \\\phi \in (0,1)\\ is analogous to a coefficient of variation.

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
brs_repar(mu = 0.5, phi = 0.2, repar = 2)
#>   shape1 shape2
#> 1      2      2
```
