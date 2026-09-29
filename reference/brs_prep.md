# Pre-process analyst data for beta interval regression

Turns analyst data into the cells \\\[l_i, u_i\]\\ and censoring types
\\\delta_i\\ used by
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md) and
[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md).
Scores run over \\0, 1, \ldots, K\\ with \\K =\\ `ncuts` the maximum;
shift a scale that starts at 1 (a Likert item 1–5 becomes 0–4,
`ncuts = 4`). Four input modes are recognised per row:

1.  **Score only** (`y`): the cell of the score and \\\delta\\ from the
    score, exactly as
    [`brs_check`](https://evandeilton.github.io/betaregscale/reference/brs_check.md)
    (0 \\\to\\ 1, \\K\\ \\\to\\ 2, otherwise 3; a value in \\(0, 1)\\
    \\\to\\ 0, per observation).

2.  **Score and `delta`**: the analyst's censoring type, with the cell
    of the score.

3.  **Bounds only** (`left` and/or `right`, `y` missing): an interval,
    or a one-sided censoring when one bound is `NA`.

4.  **Score and both bounds**: the analyst's interval.

Covariate columns are kept unchanged.

## Usage

``` r
brs_prep(
  data,
  y = "y",
  delta = "delta",
  left = "left",
  right = "right",
  ncuts = 100L,
  lim = 0.5,
  interval = c("mid", "right", "left")
)
```

## Arguments

- data:

  Data frame with the response columns and covariates.

- y:

  Character: name of the score column (default `"y"`).

- delta:

  Character: name of the censoring-type column (default `"delta"`);
  values in `{0, 1, 2, 3}` or `NA`.

- left, right:

  Character: names of the lower and upper bound columns (defaults
  `"left"`, `"right"`), on the latent score scale.

- ncuts:

  Integer \\K\\: the maximum score (default 100); the scale is \\0, 1,
  \ldots, K\\ (\\K + 1\\ categories). Must be at least the largest `y`;
  bounds may reach the latent range given in Details.

- lim:

  Numeric in \\(0, 0.5\]\\: half-width of the score cell under
  `interval = "mid"` (default 0.5); see
  [`brs_check`](https://evandeilton.github.io/betaregscale/reference/brs_check.md).

- interval:

  Direction of the uncertainty interval, `"mid"` (default), `"right"` or
  `"left"`; see
  [`brs_check`](https://evandeilton.github.io/betaregscale/reference/brs_check.md).

## Value

A data frame with columns `left`, `right` (cell on \\(0, 1)\\), `yt`
(cell centre), `y` (score, or the filled latent score for rows without
one) and `delta`, followed by the covariates. Attributes `"is_prepared"`
(`TRUE`), `"ncuts"`, `"lim"` and `"interval"` are reused by
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md) and
[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md);
an explicit different value there is ignored with a warning.

## Details

A non-missing `delta` always wins. Otherwise the type comes from the
pattern of `left`, `right` and `y` (here \\K = 10\\,
`interval = "mid"`):

|              |         |      |                                       |            |
|--------------|---------|------|---------------------------------------|------------|
| `left`       | `right` | `y`  | cell on \\(0, 1)\\                    | \\\delta\\ |
| `NA`         | 3       | `NA` | \\\[\epsilon, 0.3\]\\ (below 3)       | 1          |
| 7            | `NA`    | `NA` | \\\[0.7, 1 - \epsilon\]\\ (above 7)   | 2          |
| 2            | 5       | `NA` | \\\[0.2, 0.5\]\\                      | 3          |
| 4            | 6       | 5    | \\\[0.4, 0.6\]\\ (analyst interval)   | 3          |
| -0.5         | 4       | `NA` | \\\[\epsilon, 0.4\]\\ (reaches 0)     | 1          |
| 3            | 10.5    | `NA` | \\\[0.3, 1 - \epsilon\]\\ (reaches 1) | 2          |
| (no columns) |         | 5    | \\\[0.45, 0.55\]\\ (Mode 1)           | 3          |
| `NA`         | `NA`    | 5    | 0.5 (exact reading)                   | 0          |
| `NA`         | `NA`    | 0    | \\\[\epsilon, 0.05\]\\                | 1          |
| `NA`         | `NA`    | 10   | \\\[0.95, 1 - \epsilon\]\\            | 2          |

When the data have `left`/`right` columns, a row that gives only an
interior score (both bounds `NA`) is an exact value at the cell centre,
with a density contribution; without those columns the same score is
interval-censored (Mode 1). An analyst interval that reaches 0 (or 1) on
the unit scale is left- (or right-) censored; one that reaches both
borders covers the whole scale, stays \\\delta = 3\\ and gives a
warning, since such a row carries no information about the parameters. A
row with neither a score nor a bound is an error (a `delta` alone
defines no interval).

Score-based rows (Modes 1 and 2) use the cells of `interval`, as in
[`brs_check`](https://evandeilton.github.io/betaregscale/reference/brs_check.md):
\\\[s - \mathrm{lim}, s + \mathrm{lim}\]/K\\ for `"mid"`, \\\[s, s +
1\]/(K + 1)\\ for `"right"`/`"left"`; a forced \\\delta = 1\\ keeps
\\u_s\\ and a forced \\\delta = 2\\ keeps \\l_s\\. Analyst bounds \\L\\
(Modes 3 and 4) are values of the latent score of the chosen direction
(the scale of `predict(type = "score")`) and map to \\L/K\\ (`"mid"`),
\\L/(K + 1)\\ (`"right"`) or \\(L + 1)/(K + 1)\\ (`"left"`), so that the
dissertation's intervals \\\[s - 0.5, s + 0.5\]\\, \\\[s, s + 1\]\\ and
\\\[s - 1, s\]\\ all give the cell of score \\s\\. They must lie in the
latent range of the direction: \\\[-0.5, K + 0.5\]\\ (`"mid"`), \\\[0,
K + 1\]\\ (`"right"`) or \\\[-1, K\]\\ (`"left"`); a bound outside it is
an error. Rows with no score get `y` = the latent score of their cell
centre, so that
[`model.frame()`](https://rdrr.io/r/stats/model.frame.html) keeps them.

Endpoints are clamped to \\\[\epsilon, 1 - \epsilon\]\\, \\\epsilon =
10^{-5}\\ (why this replaces the edge-effect transformation of Lopes,
2023: section 'Scale change and borders' of
[`brs_check`](https://evandeilton.github.io/betaregscale/reference/brs_check.md)).
A \\\delta = 3\\ row with `left >= right` is an error (zero probability;
use \\\delta = 0\\), and so is a cell the clamp squeezes to zero width
(the message tells whether `ncuts` is too large for the cells or the
analyst interval lies inside the clamp). As in
[`brs_check`](https://evandeilton.github.io/betaregscale/reference/brs_check.md),
values in \\(0, 1)\\ mixed with values \\\ge 1\\ give one warning.
Unusual combinations, e.g. \\\delta = 1\\ with \\y \neq 0\\, give a
warning but are kept.

## References

Lopes, J. E. (2023). *Modelos de regressao beta para dados de escala*.
Master's dissertation, Universidade Federal do Parana, Curitiba. URI:
https://hdl.handle.net/1884/86624.

## See also

[`brs_check`](https://evandeilton.github.io/betaregscale/reference/brs_check.md)
for the cells and the automatic classification;
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md) for
fitting the model.

## Examples

``` r
# Mode 1: score only; delta from the score (0 -> 1, K -> 2, else 3)
d1 <- data.frame(y = c(0, 3, 10), x = c(1.2, 0.4, -0.3))
brs_prep(d1, ncuts = 10)
#> brs_prep: n = 3 | exact = 0, left = 1, right = 1, interval = 1
#>      left   right      yt  y delta    x
#> 1 0.00001 0.05000 0.00001  0     1  1.2
#> 2 0.25000 0.35000 0.30000  3     3  0.4
#> 3 0.95000 0.99999 0.99999 10     2 -0.3

# Mode 2: score + analyst delta (the same score read as exact and as a cell)
d2 <- data.frame(y = c(4, 4), delta = c(0, 3))
brs_prep(d2, ncuts = 10)
#> brs_prep: n = 2 | exact = 1, left = 0, right = 0, interval = 1
#>   left right  yt y delta
#> 1 0.40  0.40 0.4 4     0
#> 2 0.35  0.45 0.4 4     3

# Mode 3: only left and/or right bounds (NA pattern gives delta)
d3 <- data.frame(left = c(NA, 7, 2), right = c(3, NA, 5))
brs_prep(d3, ncuts = 10)
#> brs_prep: n = 3 | exact = 0, left = 1, right = 1, interval = 1
#>    left   right   yt   y delta
#> 1 1e-05 0.30000 0.15 1.5     1
#> 2 7e-01 0.99999 0.85 8.5     2
#> 3 2e-01 0.50000 0.35 3.5     3

# Mode 4: score with analyst bounds (used as given, divided by K)
d4 <- data.frame(y = 5, left = 4, right = 6)
brs_prep(d4, ncuts = 10)
#> brs_prep: n = 1 | exact = 0, left = 0, right = 0, interval = 1
#>   left right  yt y delta
#> 1  0.4   0.6 0.5 5     3

# An analyst interval reaching a border is one-sided censoring
brs_prep(data.frame(left = c(-0.5, 3), right = c(4, 10.5)), ncuts = 10)
#> brs_prep: n = 2 | exact = 0, left = 1, right = 1, interval = 0
#>    left   right    yt    y delta
#> 1 1e-05 0.40000 0.175 1.75     1
#> 2 3e-01 0.99999 0.675 6.75     2

# A Likert item 1-5: shift it to 0-4, so that ncuts = 4 is the maximum
lk <- data.frame(item = c(1, 2, 5, 3, 4), x = c(0.1, 0.5, 0.9, 0.3, 0.7))
brs_prep(data.frame(y = lk$item - 1, x = lk$x), ncuts = 4)
#> brs_prep: n = 5 | exact = 0, left = 1, right = 1, interval = 3
#>      left   right      yt y delta   x
#> 1 0.00001 0.12500 0.00001 0     1 0.1
#> 2 0.12500 0.37500 0.25000 1     3 0.5
#> 3 0.87500 0.99999 0.99999 4     2 0.9
#> 4 0.37500 0.62500 0.50000 2     3 0.3
#> 5 0.62500 0.87500 0.75000 3     3 0.7

# Right-direction cells [s, s + 1] / 11; the choice is stored as attributes
p <- brs_prep(d1, ncuts = 10, interval = "right")
#> brs_prep: n = 3 | exact = 0, left = 1, right = 1, interval = 1
p[, c("left", "right", "delta")]
#>        left      right delta
#> 1 0.0000100 0.09090909     1
#> 2 0.2727273 0.36363636     3
#> 3 0.9090909 0.99999000     2
attributes(p)[c("ncuts", "lim", "interval")]
#> $ncuts
#> [1] 10
#> 
#> $lim
#> [1] 0.5
#> 
#> $interval
#> [1] "right"
#> 
```
