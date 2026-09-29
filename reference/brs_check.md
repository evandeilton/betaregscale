# Transform and validate a scale-derived response variable

Maps a score on \\\\0, 1, \ldots, K\\\\ (\\K =\\ `ncuts`) to a cell
\\\[l_s, u_s\]\\ of \\(0, 1)\\ and to a censoring type \\\delta\\ of the
complete likelihood (dissertation, eq. `eqn_verossimilhanca_geral`):
\\\delta = 0\\ density \\f(y)\\, \\\delta = 1\\ \\F(u)\\, \\\delta = 2\\
\\1 - F(l)\\, \\\delta = 3\\ \\F(u) - F(l)\\. A value in \\(0, 1)\\ is
exact (\\\delta = 0\\), observation by observation.

## Usage

``` r
brs_check(
  y,
  ncuts = 100L,
  lim = 0.5,
  delta = NULL,
  interval = c("mid", "right", "left")
)
```

## Arguments

- y:

  Numeric vector: the raw response. Can be either integer scores on the
  scale \\\\0, 1, \ldots, K\\\\ or continuous values already in \\(0,
  1)\\.

- ncuts:

  Integer: number of scale categories \\K\\ (default 100). Must be
  \\\geq \max(y)\\.

- lim:

  Numeric in \\(0, 0.5\]\\: half-width of the cell under
  `interval = "mid"` (default 0.5, adjacent cells touch). Values below
  0.5 give a partial coarsening (warning); ignored, with a warning, for
  `"right"`/`"left"`.

- delta:

  Integer vector or `NULL`. If `NULL` (default), censoring types are
  derived from the scores. If provided, must have the same length as `y`
  with elements in `{0, 1, 2, 3}`; it overrides the type per observation
  (see Details).

- interval:

  Direction of the uncertainty interval: `"mid"` (default), `"right"` or
  `"left"`; see the section 'Interval direction'.

## Value

A numeric matrix with \\n\\ rows and 5 columns:

- `left`:

  Lower endpoint \\l_i\\ on \\(0, 1)\\, clamped to \\\[\epsilon, 1 -
  \epsilon\]\\.

- `right`:

  Upper endpoint \\u_i\\ on \\(0, 1)\\, clamped to \\\[\epsilon, 1 -
  \epsilon\]\\.

- `yt`:

  Midpoint approximation \\y_t\\ for starting-value computation. Also
  enters the likelihood directly as the density argument for exact
  observations (\\\delta = 0\\); for censored observations only
  `left`/`right` enter the likelihood.

- `y`:

  Original response value (preserved unchanged).

- `delta`:

  Censoring indicator: 0 = exact (density), 1 = left-censored \\F(u)\\,
  2 = right-censored \\1 - F(l)\\, 3 = interval-censored \\F(u) -
  F(l)\\.

## Details

With `delta = NULL`, each value in \\(0, 1)\\ is exact and each other
value is a score with its cell and the type above; the same
per-observation rule as
[`brs_prep`](https://evandeilton.github.io/betaregscale/reference/brs_prep.md).
Input that mixes values in \\(0, 1)\\ with values \\\ge 1\\ is ambiguous
(proportions and scores side by side, or rescaled scores): a warning
says which rule applied. Half-point scores (0, 0.5, 1, ...) are scores
on a finer grid: use `y * 2` and `ncuts * 2`. A user-supplied `delta`
(the mechanism
[`brs_sim`](https://evandeilton.github.io/betaregscale/reference/brs_sim.md)
uses in Monte Carlo studies) forces the type per observation and keeps
the cell endpoints of the score:

|            |                                           |                  |
|------------|-------------------------------------------|------------------|
| \\\delta\\ | \\l_i\\                                   | \\u_i\\          |
| 0          | cell centre (or \\y\\ when in \\(0, 1)\\) | same             |
| 1          | \\\epsilon\\                              | \\u_s\\          |
| 2          | \\l_s\\                                   | \\1 - \epsilon\\ |
| 3          | \\l_s\\                                   | \\u_s\\          |

Under `"mid"` with `lim = 0.5` this is \\u_0 = 0.5 / K\\, \\l_K = (K -
0.5) / K\\ and \\\[l_s, u_s\] = \[(s - 0.5) / K, (s + 0.5) / K\]\\.
Scores outside \\\[0, K\]\\ are an error (as in
[`brs_prep`](https://evandeilton.github.io/betaregscale/reference/brs_prep.md)),
and so is a \\\delta = 3\\ observation with \\l_i = u_i\\
(zero-probability interval).

All endpoints are clamped to \\\[\epsilon, 1 - \epsilon\]\\, \\\epsilon
= 10^{-5}\\. `yt` is the cell centre (\\s / K\\ under `"mid"`, \\(s +
0.5) / (K + 1)\\ otherwise; \\y\\ itself for exact values): it is the
density argument for \\\delta = 0\\ and a point summary elsewhere;
censored contributions use only `left`/`right`.

**Interaction with the fitting pipeline**:

This function is called internally by `.extract_response()` when the
data does *not* carry the `"is_prepared"` attribute. If data has already
been processed by
[`brs_prep`](https://evandeilton.github.io/betaregscale/reference/brs_prep.md)
or by simulation with forced delta
([`brs_sim`](https://evandeilton.github.io/betaregscale/reference/brs_sim.md)
with `delta != NULL`), the pre-computed columns are used directly and
`brs_check()` is skipped.

## Interval direction

`interval` is the direction of the uncertainty interval around the score
(dissertation, "Mapeamento de intervalos para beta": \\m = \[s - 0.5,
s + 0.5\]\\, \\r = \[s, s + 1\]\\, \\l = \[s - 1, s\]\\):

|  |  |  |
|----|----|----|
| `interval` | cell of score \\s\\ | latent score |
| `"mid"` | \\\[s - \mathrm{lim}, s + \mathrm{lim}\] / K\\ | \\K y^\*\\ |
| `"right"` | \\\[s, s + 1\] / (K + 1)\\ | \\(K + 1) y^\*\\ |
| `"left"` | \\\[s, s + 1\] / (K + 1)\\ | \\(K + 1) y^\* - 1\\ |

The \\K + 1\\ cells of `"right"` and `"left"` are equal and partition
\\\[0, 1\]\\. This normalisation is a package choice that differs from
the dissertation, which divides \\r\\ and \\l\\ by \\K\\ (there \\r\\
and \\l\\ differ by \\1/K\\, and chapter 4 reports opposite intercept
biases for them); here `"right"` and `"left"` give the same likelihood
and coefficients, and differ only in how a fitted value is read back on
the score scale (one unit), so that opposite-bias signature disappears
by construction. The three modes are different coarsening models of the
same scores: their log-likelihoods are not comparable and
[`anova()`](https://rdrr.io/r/stats/anova.html) refuses to compare them.
`lim` applies to `"mid"` only.

The censoring type comes from the score, before any clamping: \\s = 0
\to \delta = 1\\ with \\u = u_0\\, \\s = K \to \delta = 2\\ with \\l =
l_K\\, otherwise \\\delta = 3\\.

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

[`brs_prep`](https://evandeilton.github.io/betaregscale/reference/brs_prep.md)
for the analyst-facing pre-processing function;
[`brs_sim`](https://evandeilton.github.io/betaregscale/reference/brs_sim.md)
for simulation with forced delta.

## Examples

``` r
# Scale data with boundary observations
y <- c(0, 3, 5, 7, 9, 10)
brs_check(y, ncuts = 10)
#>         left   right      yt  y delta
#> [1,] 0.00001 0.05000 0.00001  0     1
#> [2,] 0.25000 0.35000 0.30000  3     3
#> [3,] 0.45000 0.55000 0.50000  5     3
#> [4,] 0.65000 0.75000 0.70000  7     3
#> [5,] 0.85000 0.95000 0.90000  9     3
#> [6,] 0.95000 0.99999 0.99999 10     2

# Right-direction intervals: cells [s, s + 1] / 11
brs_check(y, ncuts = 10, interval = "right")
#>           left      right         yt  y delta
#> [1,] 0.0000100 0.09090909 0.04545455  0     1
#> [2,] 0.2727273 0.36363636 0.31818182  3     3
#> [3,] 0.4545455 0.54545455 0.50000000  5     3
#> [4,] 0.6363636 0.72727273 0.68181818  7     3
#> [5,] 0.8181818 0.90909091 0.86363636  9     3
#> [6,] 0.9090909 0.99999000 0.95454545 10     2

# Force all observations to be exact (delta = 0)
brs_check(y, ncuts = 10, delta = rep(0L, length(y)))
#>         left   right      yt  y delta
#> [1,] 0.00001 0.00001 0.00001  0     0
#> [2,] 0.30000 0.30000 0.30000  3     0
#> [3,] 0.50000 0.50000 0.50000  5     0
#> [4,] 0.70000 0.70000 0.70000  7     0
#> [5,] 0.90000 0.90000 0.90000  9     0
#> [6,] 0.99999 0.99999 0.99999 10     0

# Force delta = 1 on non-boundary observations: u = (y + 0.5) / K
y2 <- c(30, 60)
brs_check(y2, ncuts = 100, delta = c(1L, 1L))
#>       left right  yt  y delta
#> [1,] 1e-05 0.305 0.3 30     1
#> [2,] 1e-05 0.605 0.6 60     1
```
