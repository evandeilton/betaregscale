# Map scale scores to intervals on (0, 1) and censoring types

Maps a score on \\\\0, 1, \ldots, K\\\\ (\\K =\\ `ncuts`, the maximum
score; \\K + 1\\ categories) to a cell \\\[l_s, u_s\]\\ of \\(0, 1)\\
and to a censoring type \\\delta\\ of the complete likelihood (Lopes,
2023, eq. `eqn_verossimilhanca_geral`; see
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md)):
\\\delta = 0\\ density \\f(y)\\, \\\delta = 1\\ \\F(u)\\, \\\delta = 2\\
\\1 - F(l)\\, \\\delta = 3\\ \\F(u) - F(l)\\. A value in \\(0, 1)\\ is
exact (\\\delta = 0\\), observation by observation. A scale that starts
at 1 (e.g. a Likert item 1–5) must be shifted to start at 0 (0–4,
`ncuts = 4`): its lowest category is then the left-censored border.

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

  Numeric vector: scores on \\\\0, 1, \ldots, K\\\\ or values in \\(0,
  1)\\.

- ncuts:

  Integer \\K\\: the maximum score (default 100); the scale is \\0, 1,
  \ldots, K\\, with \\K + 1\\ categories. Must be \\\geq \max(y)\\.

- lim:

  Numeric in \\(0, 0.5\]\\: half-width of the cell under
  `interval = "mid"` (default 0.5, adjacent cells touch). Values below
  0.5 give a partial coarsening (warning); ignored, with a warning, for
  `"right"`/`"left"`.

- delta:

  Integer vector or `NULL`. If `NULL` (default), censoring types are
  derived from the scores. If provided, it must have the length of `y`
  with elements in `{0, 1, 2, 3}`; it overrides the type per observation
  (see Details).

- interval:

  Direction of the uncertainty interval: `"mid"` (default), `"right"` or
  `"left"`; see the section 'Interval direction'.

## Value

A numeric matrix with one row per observation and columns `left`
(\\l_i\\), `right` (\\u_i\\), `yt` (cell centre), `y` (the input) and
`delta`.

## Details

With `delta = NULL`, each value in \\(0, 1)\\ is exact and each other
value is a score with its cell and the type above; this is the rule of
[`brs_prep`](https://evandeilton.github.io/betaregscale/reference/brs_prep.md)
too. Input that mixes values in \\(0, 1)\\ with values \\\ge 1\\ is
ambiguous (proportions next to scores, or rescaled scores) and gives one
warning; rescale it if the values in \\(0, 1)\\ are meant as scores
(half-point scores 0, 0.5, 1, ...: use `y * 2` and `ncuts * 2`). A
user-supplied `delta` (the mechanism
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

A forced \\\delta = 1\\ or 2 uses the cell of the observed score as the
censoring threshold, so the threshold depends on \\Y\\: applied to every
observation this is informative censoring and the likelihood has no
finite maximum
([`brs_sim`](https://evandeilton.github.io/betaregscale/reference/brs_sim.md)
warns). Under `"mid"` with `lim = 0.5`: \\u_0 = 0.5/K\\, \\l_K = (K -
0.5)/K\\ and \\\[l_s, u_s\] = \[(s - 0.5)/K, (s + 0.5)/K\]\\. Scores
outside \\\[0, K\]\\ are an error (as in
[`brs_prep`](https://evandeilton.github.io/betaregscale/reference/brs_prep.md)),
and so is a \\\delta = 3\\ observation with \\l_i \ge u_i\\
(zero-probability interval). When the clamp squeezes a cell to zero
width (a cell narrower than \\10^{-5}\\ at a border, i.e. `ncuts` too
large or `lim` too small), the error says so.

`yt` is the cell centre: \\s/K\\ under `"mid"` (so the border scores sit
at \\\epsilon\\ and \\1 - \epsilon\\) and \\(s + 0.5)/(K + 1)\\
otherwise; exact values are kept. It is the density argument for
\\\delta = 0\\ and the point used by the midpoint residuals; censored
contributions do not use it.

Data carrying the `"is_prepared"` attribute (from
[`brs_prep`](https://evandeilton.github.io/betaregscale/reference/brs_prep.md)
or
[`brs_sim`](https://evandeilton.github.io/betaregscale/reference/brs_sim.md))
are used by
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md) as
they are; `brs_check()` runs only for raw scores.

## Interval direction

`interval` is the direction of the uncertainty interval around the score
(Lopes, 2023, "Mapeamento de intervalos para beta": \\m = \[s - 0.5, s +
0.5\]\\, \\r = \[s, s + 1\]\\, \\l = \[s - 1, s\]\\).

|  |  |  |
|----|----|----|
| `interval` | cell of score \\s\\ | latent score of \\y^\* \in (0, 1)\\ |
| `"mid"` | \\\[s - \mathrm{lim}, s + \mathrm{lim}\]/K\\ | \\K y^\*\\ |
| `"right"` | \\\[s, s + 1\]/(K + 1)\\ | \\(K + 1) y^\*\\ |
| `"left"` | \\\[s, s + 1\]/(K + 1)\\ | \\(K + 1) y^\* - 1\\ |

The latent score is the back-transformation used by
`predict(type = "score")`: it lies in \\\[s - 0.5, s + 0.5\]\\, \\\[s,
s + 1\]\\ or \\\[s - 1, s\]\\ when \\y^\*\\ is in the cell of \\s\\. The
\\K + 1\\ cells of `"right"` and `"left"` are equal and partition \\\[0,
1\]\\. This is a package choice that differs from the dissertation,
which divides \\r\\ and \\l\\ by \\K\\; there the two directions differ
by \\1/K\\ and chapter 4 reports opposite intercept biases for them,
while here they give the same likelihood and the same coefficients. The
three directions are different coarsening models of the same scores:
their log-likelihoods are not comparable, and
[`anova.brs`](https://evandeilton.github.io/betaregscale/reference/anova.brs.md)
refuses such comparisons. `lim` applies to `"mid"` only.

## Scale change and borders

Lopes (2023, "Mudanca de escala intervalar") rescales scores by the
range and handles the edge effect either with the transformation \\y^\*
= \\y(n - 1)/R + 1/2\\/n\\ (Smithson and Verkuilen, 2006), which depends
on the sample size \\n\\ and the observed range \\R\\, or by moving 0
and 1 inwards by \\\zeta = 10^{-4}\\. The package does neither: each
score is coarsened into its cell, using the scale maximum \\K\\ rather
than the observed range; the border scores 0 and \\K\\ become left- and
right-censored observations with contributions \\F(u_0)\\ and \\1 -
F(l_K)\\; and the endpoints are clamped to \\\[\epsilon, 1 -
\epsilon\]\\, \\\epsilon = 10^{-5}\\, only as a numerical guard. The
border cells therefore keep their full probability, no data are moved
towards 1/2, and the fit does not depend on \\n\\. The clamp does not
enter the border contributions (\\F(u_0)\\ ignores \\l_0\\, \\1 -
F(l_K)\\ ignores \\u_K\\); it matters only for exact values and for
extreme user-supplied intervals.

The censoring type is derived from the score, before any clamping: \\s =
0 \to \delta = 1\\, \\s = K \to \delta = 2\\, otherwise \\\delta = 3\\.

## References

Lopes, J. E. (2023). *Modelos de regressao beta para dados de escala*.
Master's dissertation, Universidade Federal do Parana, Curitiba. URI:
https://hdl.handle.net/1884/86624.

Smithson, M., and Verkuilen, J. (2006). A better lemon squeezer?
Maximum-likelihood regression with beta-distributed dependent variables.
*Psychological Methods*, **11**(1), 54–71.
[doi:10.1037/1082-989X.11.1.54](https://doi.org/10.1037/1082-989X.11.1.54)

Hjermstad, M. J., Fayers, P. M., Haugen, D. F., et al. (2011). Studies
comparing Numerical Rating Scales, Verbal Rating Scales, and Visual
Analogue Scales for assessment of pain intensity in adults: a systematic
literature review. Journal of Pain and Symptom Management, 41(6),
1073-1093.
[doi:10.1016/j.jpainsymman.2010.08.016](https://doi.org/10.1016/j.jpainsymman.2010.08.016)

## See also

[`brs_prep`](https://evandeilton.github.io/betaregscale/reference/brs_prep.md)
for analyst-supplied censoring and intervals;
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md);
[`brs_sim`](https://evandeilton.github.io/betaregscale/reference/brs_sim.md).

## Examples

``` r
# NRS-11 scores (K = 10) under the three interval directions
s <- 0:10
brs_check(s, ncuts = 10)                      # mid:   [s - 0.5, s + 0.5] / 10
#>          left   right      yt  y delta
#>  [1,] 0.00001 0.05000 0.00001  0     1
#>  [2,] 0.05000 0.15000 0.10000  1     3
#>  [3,] 0.15000 0.25000 0.20000  2     3
#>  [4,] 0.25000 0.35000 0.30000  3     3
#>  [5,] 0.35000 0.45000 0.40000  4     3
#>  [6,] 0.45000 0.55000 0.50000  5     3
#>  [7,] 0.55000 0.65000 0.60000  6     3
#>  [8,] 0.65000 0.75000 0.70000  7     3
#>  [9,] 0.75000 0.85000 0.80000  8     3
#> [10,] 0.85000 0.95000 0.90000  9     3
#> [11,] 0.95000 0.99999 0.99999 10     2
brs_check(s, ncuts = 10, interval = "right")  # right: [s, s + 1] / 11
#>             left      right         yt  y delta
#>  [1,] 0.00001000 0.09090909 0.04545455  0     1
#>  [2,] 0.09090909 0.18181818 0.13636364  1     3
#>  [3,] 0.18181818 0.27272727 0.22727273  2     3
#>  [4,] 0.27272727 0.36363636 0.31818182  3     3
#>  [5,] 0.36363636 0.45454545 0.40909091  4     3
#>  [6,] 0.45454545 0.54545455 0.50000000  5     3
#>  [7,] 0.54545455 0.63636364 0.59090909  6     3
#>  [8,] 0.63636364 0.72727273 0.68181818  7     3
#>  [9,] 0.72727273 0.81818182 0.77272727  8     3
#> [10,] 0.81818182 0.90909091 0.86363636  9     3
#> [11,] 0.90909091 0.99999000 0.95454545 10     2
# "left" has the same cells as "right"; only the latent score read back differs
identical(brs_check(s, 10, interval = "left"), brs_check(s, 10, interval = "right"))
#> [1] TRUE

# Borders: score 0 -> delta 1 (F(u)), score 10 -> delta 2 (1 - F(l))
brs_check(c(0, 10), ncuts = 10)[, c("left", "right", "delta")]
#>         left   right delta
#> [1,] 0.00001 0.05000     1
#> [2,] 0.95000 0.99999     2

# A Likert item 1-5: shift to 0-4 so that K = ncuts = 4 is its maximum
likert <- c(1, 2, 3, 5, 4)
brs_check(likert - 1, ncuts = 4)[, c("left", "right", "delta")]
#>         left   right delta
#> [1,] 0.00001 0.12500     1
#> [2,] 0.12500 0.37500     3
#> [3,] 0.37500 0.62500     3
#> [4,] 0.87500 0.99999     2
#> [5,] 0.62500 0.87500     3
# Unshifted, category 1 would be an interior cell and 0 a category never used
brs_check(likert, ncuts = 5)[, "delta"]
#> [1] 3 3 3 2 3

# Values already in (0, 1) are exact observations (delta 0)
brs_check(c(0.12, 0.5, 0.97))
#> Response is already on the unit interval (0, 1); treating as uncensored (exact) observations.
#>      left right   yt    y delta
#> [1,] 0.12  0.12 0.12 0.12     0
#> [2,] 0.50  0.50 0.50 0.50     0
#> [3,] 0.97  0.97 0.97 0.97     0
# Per observation: 0.3 is exact, 5 and 10 are scores (with a warning)
brs_check(c(0.3, 5, 10), ncuts = 10)
#> Warning: The response mixes values in (0, 1) with values >= 1: values in (0, 1) are taken as exact (delta = 0), the others as scores on 0..ncuts. Rescale the data if the (0, 1) values are scores (half-point scores: use y * 2 and ncuts * 2).
#>      left   right      yt    y delta
#> [1,] 0.30 0.30000 0.30000  0.3     0
#> [2,] 0.45 0.55000 0.50000  5.0     3
#> [3,] 0.95 0.99999 0.99999 10.0     2

# A forced delta keeps the cell of the score
brs_check(c(30, 60), ncuts = 100, delta = c(1L, 2L))
#>         left   right  yt  y delta
#> [1,] 0.00001 0.30500 0.3 30     1
#> [2,] 0.59500 0.99999 0.6 60     2
```
