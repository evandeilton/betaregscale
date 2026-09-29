# Fit a beta interval regression model

Unified interface that dispatches to
[`brs_fit_fixed`](https://evandeilton.github.io/betaregscale/reference/brs_fit_fixed.md)
(fixed dispersion) or
[`brs_fit_var`](https://evandeilton.github.io/betaregscale/reference/brs_fit_var.md)
(variable dispersion) based on the formula structure.

## Usage

``` r
brs(
  formula,
  data,
  link = NULL,
  link_phi = NULL,
  ncuts = NULL,
  lim = NULL,
  repar = 2L,
  method = c("BFGS", "L-BFGS-B"),
  hessian_method = c("numDeriv", "optim"),
  interval = NULL
)
```

## Arguments

- formula:

  A [`Formula`](https://rdrr.io/pkg/Formula/man/Formula.html)-style
  formula with two parts: `y ~ x1 + x2 | z1 + z2`.

- data:

  Data frame.

- link:

  Link for the first parameter (the mean under `repar = 1, 2`; the shape
  \\p\\ under `repar = 0`). `NULL` (default) selects the link implied by
  `repar`: `"logit"` for `repar = 1, 2`, `"log"` for `repar = 0`. See
  the 'Reparameterizations and links' section of `brs` for the
  admissible values.

- link_phi:

  Link for the second parameter. `NULL` (default) selects `"logit"` for
  `repar = 2` (dispersion on \\(0, 1)\\) and `"log"` for `repar = 0, 1`
  (positive shape/precision).

- ncuts:

  Number of scale categories. `NULL` (default) uses the value stored by
  [`brs_prep`](https://evandeilton.github.io/betaregscale/reference/brs_prep.md)
  in `attr(data, "ncuts")`, or 100 when `data` was not prepared. A value
  that differs from the stored one is ignored with a warning (the
  endpoints were built with the stored value).

- lim:

  Half-width of the score cell in \\(0, 0.5\]\\ (`interval = "mid"`
  only). `NULL` (default) uses `attr(data, "lim")` from
  [`brs_prep`](https://evandeilton.github.io/betaregscale/reference/brs_prep.md),
  or 0.5; same rule as `ncuts`. Values below 0.5 warn (partial
  coarsening).

- repar:

  Reparameterization scheme (default 2); see
  [`brs_repar`](https://evandeilton.github.io/betaregscale/reference/brs_repar.md).

- method:

  Optimization method (default `"BFGS"`).

- hessian_method:

  Character: `"numDeriv"` or `"optim"`.

- interval:

  Direction of the uncertainty interval, `"mid"`, `"right"` or `"left"`
  (see
  [`brs_check`](https://evandeilton.github.io/betaregscale/reference/brs_check.md)).
  `NULL` (default) uses `attr(data, "interval")` from
  [`brs_prep`](https://evandeilton.github.io/betaregscale/reference/brs_prep.md),
  or `"mid"`; same rule as `ncuts`.

## Value

An object of class `"brs"`.

## Details

If the formula contains a `|` separator (e.g., `y ~ x1 + x2 | z1`), the
variable-dispersion model is fitted; otherwise, a fixed-dispersion model
is used.

## Reparameterizations and links

The three schemes of
[`brs_repar`](https://evandeilton.github.io/betaregscale/reference/brs_repar.md)
model different parameters, so the admissible links differ: parameters
on \\(0, 1)\\ use a `(0, 1)`-link, parameters on \\(0, \infty)\\ use
`"log"` or `"sqrt"`. `link = NULL` and `link_phi = NULL` (the defaults)
select the first entry of each cell; any other combination is rejected
with an error.

|  |  |  |
|----|----|----|
| `repar` | `link` (first parameter) | `link_phi` (second parameter) |
| 0 (shapes \\p, q\\) | `log`, `sqrt` | `log`, `sqrt` |
| 1 (mean, precision) | `logit`, `probit`, `cauchit`, `cloglog` | `log`, `sqrt` |
| 2 (mean, dispersion) | `logit`, `probit`, `cauchit`, `cloglog` | `logit`, `probit`, `cauchit`, `cloglog` |

`"identity"`, `"inverse"` and `"1/mu^2"` are not accepted for positive
parameters (their inverse does not map the real line onto \\(0,
\infty)\\). With `"sqrt"` the inverse link is flat for \\\eta \le 0\\; a
warning is issued after the fit when a fitted linear predictor lies on
that plateau.

Under `repar = 0` the fitted object stores the shape \\p\\ in `hatmu`
(and `predict(type = "link")` is its linear predictor), while
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
`predict(type = "response")`, residuals and marginal effects use the
mean \\E\[Y\] = p / (p + q)\\.

## Interval direction

`interval` selects how a score \\s\\ is coarsened into a cell of \\(0,
1)\\: `"mid"` \\\[s - \mathrm{lim}, s + \mathrm{lim}\] / K\\ (default),
`"right"` and `"left"` \\\[s, s + 1\] / (K + 1)\\ (the dissertation's
\\r\\ and \\l\\ directions; equal cells, a package normalisation).
`"right"` and `"left"` give the same likelihood and coefficients and
differ only in the latent score read back by `predict(type = "score")`
(one unit). The modes are different coarsening models, so their
log-likelihoods are not comparable; see
[`brs_check`](https://evandeilton.github.io/betaregscale/reference/brs_check.md)
and
[`anova.brs`](https://evandeilton.github.io/betaregscale/reference/anova.brs.md).

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
# Fixed dispersion
fit1 <- brs(y ~ x1, data = prep)
print(fit1)
#> 
#> Call:
#> brs(formula = y ~ x1, data = prep)
#> 
#> Coefficients (mean model with logit link):
#> (Intercept)          x1 
#>      0.2551     -0.2202 
#> 
#> Phi coefficients (precision model with logit link):
#>   (phi) 
#> -0.3929 
#> 
# Variable dispersion
fit2 <- brs(y ~ x1 | x2, data = prep)
print(fit2)
#> 
#> Call:
#> brs(formula = y ~ x1 | x2, data = prep)
#> 
#> Coefficients (mean model with logit link):
#> (Intercept)          x1 
#>      0.2732     -0.2310 
#> 
#> Phi coefficients (precision model with logit link):
#> (Intercept)          x2 
#>     -0.3789     -0.0288 
#> 
# }
```
