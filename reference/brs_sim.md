# Simulate data from beta interval models

Simulates scores on \\0, 1, \ldots, K\\ (\\K =\\ `ncuts`) from a fixed-
or variable-dispersion beta regression, coarsened exactly as the
likelihood of
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md)
assumes. The output can be passed to
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md) as
it is.

## Usage

``` r
brs_sim(
  formula,
  data,
  beta,
  phi = 1/5,
  zeta = NULL,
  link = NULL,
  link_phi = NULL,
  ncuts = 100L,
  lim = 0.5,
  repar = 2L,
  delta = NULL,
  interval = c("mid", "right", "left")
)
```

## Arguments

- formula:

  Model formula with one (mean) or two parts (mean `|` precision).

- data:

  Data frame with the covariates.

- beta:

  Coefficients of the first parameter, on its link scale.

- phi:

  Second parameter on its link scale (one-part formulas).

- zeta:

  Coefficients of the second parameter, on its link scale (two-part
  formulas).

- link, link_phi:

  Links; `NULL` (default) selects those implied by `repar` (see
  [`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md)).

- ncuts:

  Integer \\K\\: the maximum score (default 100); scores run over \\0,
  \ldots, K\\.

- lim:

  Half-width of the cell under `"mid"` (default 0.5).

- repar:

  Parameterisation (0, 1 or 2); see
  [`brs_repar`](https://evandeilton.github.io/betaregscale/reference/brs_repar.md).

- delta:

  `NULL` (default) or one censoring type (0, 1, 2 or 3) forced on every
  row; see Details.

- interval:

  `"mid"` (default), `"right"` or `"left"`.

## Value

A data frame with columns `left`, `right`, `yt`, `y` (score), `delta`
and the non-intercept columns of the model matrices, with the attributes
`"is_prepared"`, `"ncuts"`, `"lim"` and `"interval"` that
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md)
reuses.

## Details

`formula` follows
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md): a
one-part formula (`~ x1 + x2`) uses the scalar `phi`, a two-part formula
(`~ x1 + x2 | z1`) the coefficients `zeta`. A left-hand side is ignored.
For each row, \\Y \sim \mathrm{Beta}(a, b)\\ with the shapes of
[`brs_repar`](https://evandeilton.github.io/betaregscale/reference/brs_repar.md),
and the score is \\s = \mathrm{round}(K Y)\\ under `"mid"` (cells
centred on the score; exact only for `lim = 0.5`) or \\s = \lfloor
(K + 1) Y \rfloor\\ under `"right"`/`"left"` (\\K + 1\\ equal cells).
The cells and censoring types then come from
[`brs_check`](https://evandeilton.github.io/betaregscale/reference/brs_check.md)
with the same `interval`: 0 is left-censored, \\K\\ right-censored, the
rest interval-censored. Censoring is therefore non-informative: the
cells are fixed and do not depend on \\Y\\.

`delta` overrides the type for every row: `0` keeps the continuous \\Y\\
as exact values; `3` keeps the scores away from the borders (\\1 \le s
\le K - 1\\); `1` or `2` censors every row on one side at the cell of
its own score. The last two make the threshold depend on \\Y\\
(informative censoring): the likelihood has no finite maximum and
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md)
returns diverging coefficients, with clamp warnings. `brs_sim()` warns
in that case, and also whenever all simulated rows fall on the same
border.

A Monte Carlo study in the design of Lopes (2023, ch. 3) is a loop of
`brs_sim()` and
[`brs()`](https://evandeilton.github.io/betaregscale/reference/brs.md)
over a fixed design; see the vignette
[`vignette("brs-advanced-workflows")`](https://evandeilton.github.io/betaregscale/articles/brs-advanced-workflows.md).

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

## See also

[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md),
[`brs_check`](https://evandeilton.github.io/betaregscale/reference/brs_check.md),
[`brs_repar`](https://evandeilton.github.io/betaregscale/reference/brs_repar.md),
[`brs_bootstrap`](https://evandeilton.github.io/betaregscale/reference/brs_bootstrap.md)

## Examples

``` r
# One design, three parameterisations with their default links
set.seed(42)
d <- data.frame(x = rnorm(300))

# repar = 2 (logit/logit): mean and dispersion on (0, 1)
s2 <- brs_sim(~ x, data = d, beta = c(0.3, 0.6), phi = qlogis(0.2), ncuts = 10)
f2 <- brs(y ~ x, data = s2)
round(cbind(true = c(0.3, 0.6, qlogis(0.2)), est = coef(f2)), 2)
#>              true   est
#> (Intercept)  0.30  0.34
#> x            0.60  0.59
#> (phi)       -1.39 -1.42

# repar = 1 (logit/log): precision 20
s1 <- brs_sim(~ x, data = d, beta = c(0.3, 0.6), phi = log(20), ncuts = 10,
              repar = 1)
f1 <- brs(y ~ x, data = s1, repar = 1)
round(cbind(true = c(0.3, 0.6, log(20)), est = coef(f1)), 2)
#>             true  est
#> (Intercept)  0.3 0.28
#> x            0.6 0.59
#> (phi)        3.0 2.99

# repar = 0 (log/log): shapes p = 3 exp(0.3 x), q = 2
s0 <- brs_sim(~ x, data = d, beta = c(log(3), 0.3), phi = log(2), ncuts = 10,
              repar = 0)
f0 <- brs(y ~ x, data = s0, repar = 0)
round(cbind(true = c(log(3), 0.3, log(2)), est = coef(f0)), 2)
#>             true  est
#> (Intercept) 1.10 1.05
#> x           0.30 0.28
#> (phi)       0.69 0.63

# The output is ready for brs(): endpoints, delta and attributes
head(s2)
#>   left right  yt y delta          x
#> 1 0.75  0.85 0.8 8     3  1.3709584
#> 2 0.65  0.75 0.7 7     3 -0.5646982
#> 3 0.55  0.65 0.6 6     3  0.3631284
#> 4 0.35  0.45 0.4 4     3  0.6328626
#> 5 0.65  0.75 0.7 7     3  0.4042683
#> 6 0.55  0.65 0.6 6     3 -0.1061245
attributes(s2)[c("ncuts", "lim", "interval")]
#> $ncuts
#> [1] 10
#> 
#> $lim
#> [1] 0.5
#> 
#> $interval
#> [1] "mid"
#> 

# Right-direction cells: scores drawn as floor(11 * Y)
sr <- brs_sim(~ x, data = d, beta = c(0.3, 0.6), phi = qlogis(0.2),
              ncuts = 10, interval = "right")
table(sr$delta)
#> 
#>   1   2   3 
#>  15  25 260 

# Forcing delta = 1 on every row is informative censoring: a warning
sl <- brs_sim(~ x, data = d, beta = c(0.3, 0.6), phi = qlogis(0.2),
              ncuts = 10, delta = 1)
#> Warning: delta = 1 censors every observation on the same side at the cell of its own value (informative censoring): brs() has no finite MLE for these data (estimates diverge).
table(sl$delta)
#> 
#>   1 
#> 300 
```
