# Transform and validate a scale-derived response variable

Takes a discrete (or continuous) response on the scale \\0, 1, \ldots,
K\\ (where \\K =\\ `ncuts`) and converts it to a pair of interval
endpoints on the open unit interval \\(0, 1)\\. Each observation is
classified into one of four censoring types following the complete
likelihood of Lopes (2024, Eq. 2.24):

- \\\delta = 0\\:

  Uncensored (exact): the observation is a continuous value already in
  \\(0, 1)\\. The likelihood contribution is the density \\f(y_i \|
  \theta)\\.

- \\\delta = 1\\:

  Left-censored: the observation equals the scale minimum (\\y = 0\\).
  We only know that the latent value is below some upper bound. The
  contribution is \\F(u_i \| \theta)\\.

- \\\delta = 2\\:

  Right-censored: the observation equals the scale maximum (\\y = K\\).
  The contribution is \\1 - F(l_i \| \theta)\\.

- \\\delta = 3\\:

  Interval-censored: the standard case for scale data. The contribution
  is \\F(u_i \| \theta) - F(l_i \| \theta)\\.

## Usage

``` r
check_response(y, type = "m", ncuts = 100L, lim = 0.5, delta = NULL)
```

## Arguments

- y:

  Numeric vector — the raw response.

- type:

  **Deprecated.** Character: interval type. `"m"` = midpoint (default),
  `"l"` = left-aligned, `"r"` = right-aligned. Use
  [`bs_prepare`](https://evandeilton.github.io/betaregscale/reference/bs_prepare.md)
  to control interval geometry instead.

- ncuts:

  Integer: number of scale categories (default 100).

- lim:

  Numeric: half-width of the uncertainty region (default 0.5).

- delta:

  Integer vector or `NULL`. If `NULL` (default), censoring types are
  inferred automatically from boundary rules. If provided, must have the
  same length as `y` with values in `{0, 1, 2, 3}`. When provided,
  overrides the automatic classification.

## Value

A matrix with columns `left`, `right`, `yt` (midpoint approximation),
`y` (original value), and `delta` (censoring indicator: 0 = exact, 1 =
left, 2 = right, 3 = interval).

## Details

If the input is already in \\(0, 1)\\ (i.e., all values satisfy \\0 \< y
\< 1\\), observations are treated as uncensored (\\\delta = 0\\), unless
`delta` is provided.

If `delta` is supplied, the user-provided censoring indicators are used
directly, overriding the automatic boundary-based classification rules.

For scale (integer) data (when `delta = NULL`):

- \\y = 0\\: left-censored. Upper bound is \\u = \mathrm{lim} / K\\.

- \\y = K\\: right-censored. Lower bound is \\l = (K - \mathrm{lim}) /
  K\\.

- \\0 \< y \< K\\: interval-censored with endpoints determined by the
  interval `type`.

Three interval types (Section 2.7, Figure 4 of Lopes, 2024):

- `"m"` (midpoint):

  Symmetric: \\l = (y - \mathrm{lim})/K\\, \\u = (y + \mathrm{lim})/K\\.

- `"l"` (left):

  Shifted left: \\l = (y - 2\mathrm{lim})/K\\, \\u = y/K\\.

- `"r"` (right):

  Shifted right: \\l = y/K\\, \\u = (y + 2\mathrm{lim})/K\\.

All endpoints are clamped to \\\[\epsilon, 1 - \epsilon\]\\ with
\\\epsilon = 10^{-5}\\ to avoid boundary issues in the beta likelihood.

## Examples

``` r
# Scale data with boundary observations
y <- c(0, 3, 5, 7, 9, 10)
check_response(y, ncuts = 10)
#>         left   right      yt  y delta
#> [1,] 0.00001 0.05000 0.00001  0     1
#> [2,] 0.25000 0.35000 0.30000  3     3
#> [3,] 0.45000 0.55000 0.50000  5     3
#> [4,] 0.65000 0.75000 0.70000  7     3
#> [5,] 0.85000 0.95000 0.90000  9     3
#> [6,] 0.95000 0.99999 0.99999 10     2

# Force all observations to be exact (delta = 0)
check_response(y, ncuts = 10, delta = rep(0L, length(y)))
#>         left   right      yt  y delta
#> [1,] 0.00001 0.00001 0.00001  0     0
#> [2,] 0.30000 0.30000 0.30000  3     0
#> [3,] 0.50000 0.50000 0.50000  5     0
#> [4,] 0.70000 0.70000 0.70000  7     0
#> [5,] 0.90000 0.90000 0.90000  9     0
#> [6,] 0.99999 0.99999 0.99999 10     0
```
