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

  Numeric vector: the first parameter. The mean, in \\(0, 1)\\, for
  `repar = 1, 2`; the shape \\p \> 0\\ for `repar = 0`.

- phi:

  Numeric vector (or scalar): the second parameter (precision \\\phi \>
  0\\, dispersion \\\phi \in (0, 1)\\, or shape \\q \> 0\\).

- repar:

  Integer (0, 1, or 2) selecting the scheme.

## Value

A `data.frame` with columns `shape1` and `shape2`.

## Details

The three schemes are:

- `repar = 0`:

  Shapes: \\a = p,\\ b = q\\ with \\p, q \> 0\\ (the first argument is
  \\p\\, the second \\q\\). Both parameters live on \\(0, \infty)\\ and
  the mean is \\E\[Y\] = p / (p + q)\\. Regression directly on the
  shapes is a package extension: the dissertation presents the \\(p,
  q)\\ form (its eq. `eqn_beta_p1`) and builds the regression models on
  the two reparameterizations below.

- `repar = 1`:

  Ferrari–Cribari-Neto (dissertation "parametrização 1", eq.
  `eqn_beta_p2`): \\a = \mu\phi,\\ b = (1 - \mu)\phi\\, where \\\mu \in
  (0, 1)\\ is the mean and \\\phi \> 0\\ acts as a precision parameter.

- `repar = 2`:

  Mean–dispersion (dissertation "parametrização 2", eq. `eqn_beta_p3`):
  \\a = \mu(1-\phi)/\phi,\\ b = (1-\mu)(1-\phi)/\phi\\, where \\\mu \in
  (0, 1)\\ is the mean and \\\phi \in (0,1)\\ is a dispersion parameter
  (\\Var\[Y\] = \phi\\\mu(1-\mu)\\).

The admissible link functions follow from these domains; see the
'Reparameterizations and links' section of
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md).

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
brs_repar(mu = 0.5, phi = 0.2, repar = 2)
#>   shape1 shape2
#> 1      2      2
```
