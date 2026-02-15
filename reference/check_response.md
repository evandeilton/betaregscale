# Transform and validate a scale-derived response variable

Takes a discrete (or continuous) response on the scale \\0, 1, \ldots,
K\\ (where \\K =\\ `ncuts`) and converts it to a pair of interval
endpoints on the open unit interval \\(0, 1)\\. Each observation is
classified into one of four censoring types following the complete
likelihood of Lopes (2024, Eq. 2.24):

- \\\delta = 0\\:

  Uncensored (exact): the observation is a continuous value already in
  \\(0, 1)\\. The likelihood contribution is the density \\f(y_i \|
  \theta)\\. Endpoints: \\l_i = u_i = y_i\\ (or \\y_i / K\\ when on the
  scale).

- \\\delta = 1\\:

  Left-censored: the latent value is below some upper bound \\u_i\\. The
  contribution is \\F(u_i \| \theta)\\. When the observation is at the
  scale minimum (\\y = 0\\), the upper bound is \\u_i = \mathrm{lim} /
  K\\. When the user forces \\\delta = 1\\ on a non-boundary observation
  (\\y \neq 0\\), the upper bound is \\u_i = (y + \mathrm{lim}) / K\\,
  preserving observation- specific variation. In both cases \\l_i =
  \epsilon\\.

- \\\delta = 2\\:

  Right-censored: the latent value is above some lower bound \\l_i\\.
  The contribution is \\1 - F(l_i \| \theta)\\. When the observation is
  at the scale maximum (\\y = K\\), the lower bound is \\l_i = (K -
  \mathrm{lim}) / K\\. When the user forces \\\delta = 2\\ on a
  non-boundary observation (\\y \neq K\\), the lower bound is \\l_i =
  (y - \mathrm{lim}) / K\\, preserving observation- specific variation.
  In both cases \\u_i = 1 - \epsilon\\.

- \\\delta = 3\\:

  Interval-censored: the standard case for scale data. The contribution
  is \\F(u_i \| \theta) - F(l_i \| \theta)\\. Endpoints depend on the
  interval `type` (see Details).

## Usage

``` r
check_response(y, type = "m", ncuts = 100L, lim = 0.5, delta = NULL)
```

## Arguments

- y:

  Numeric vector: the raw response. Can be either integer scores on the
  scale \\\\0, 1, \ldots, K\\\\ or continuous values already in \\(0,
  1)\\.

- type:

  **Deprecated.** Character: interval type for \\\delta = 3\\
  observations. `"m"` = midpoint (default), `"l"` = left-aligned, `"r"`
  = right-aligned. Use
  [`bs_prepare`](https://evandeilton.github.io/betaregscale/reference/bs_prepare.md)
  to control interval geometry instead.

- ncuts:

  Integer: number of scale categories \\K\\ (default 100). Must be
  \\\geq \max(y)\\.

- lim:

  Numeric: half-width \\h\\ of the uncertainty region (default 0.5).
  Controls the width of the interval around each scale point.

- delta:

  Integer vector or `NULL`. If `NULL` (default), censoring types are
  inferred automatically from the boundary rules described above.

  If provided, must have the same length as `y`, with every element in
  `{0, 1, 2, 3}`. The supplied values override the automatic
  classification on a per-observation basis, and the endpoint formulas
  adapt to non-boundary observations as described in the table above.

  This parameter is used internally by the simulation functions when the
  analyst forces a specific censoring type (e.g.,
  `betaregscale_simulate(..., delta = 2)`).

## Value

A numeric matrix with \\n\\ rows and 5 columns:

- `left`:

  Lower endpoint \\l_i\\ on \\(0, 1)\\, clamped to \\\[\epsilon, 1 -
  \epsilon\]\\.

- `right`:

  Upper endpoint \\u_i\\ on \\(0, 1)\\, clamped to \\\[\epsilon, 1 -
  \epsilon\]\\.

- `yt`:

  Midpoint approximation \\y_t\\ for starting-value computation (does
  not enter the likelihood).

- `y`:

  Original response value (preserved unchanged).

- `delta`:

  Censoring indicator: 0 = exact (density), 1 = left-censored \\F(u)\\,
  2 = right-censored \\1 - F(l)\\, 3 = interval-censored \\F(u) -
  F(l)\\.

## Details

**Automatic classification** (`delta = NULL`):

If the entire input vector is already in \\(0, 1)\\ (i.e., all values
satisfy \\0 \< y \< 1\\), all observations are treated as uncensored
(\\\delta = 0\\).

Otherwise, for scale (integer) data:

- \\y = 0\\: left-censored (\\\delta = 1\\).

- \\y = K\\: right-censored (\\\delta = 2\\).

- \\0 \< y \< K\\: interval-censored (\\\delta = 3\\).

**User-supplied delta** (`delta` vector):

When the `delta` argument is provided, the user-supplied censoring
indicators override the automatic boundary-based rules on a
per-observation basis. This is the mechanism used by
[`betaregscale_simulate`](https://evandeilton.github.io/betaregscale/reference/betaregscale_simulate.md)
when the analyst forces a specific censoring type in Monte Carlo
studies.

The endpoint formulas for each delta value are:

|            |                       |                             |                             |
|------------|-----------------------|-----------------------------|-----------------------------|
| \\\delta\\ | Condition             | \\l_i\\ (left)              | \\u_i\\ (right)             |
| 0          | \\y \in (0, 1)\\      | \\y\\                       | \\y\\                       |
| 0          | \\y\\ on scale        | \\y / K\\                   | \\y / K\\                   |
| 1          | \\y = 0\\ (boundary)  | \\\epsilon\\                | \\\mathrm{lim} / K\\        |
| 1          | \\y \neq 0\\ (forced) | \\\epsilon\\                | \\(y + \mathrm{lim}) / K\\  |
| 2          | \\y = K\\ (boundary)  | \\(K - \mathrm{lim}) / K\\  | \\1 - \epsilon\\            |
| 2          | \\y \neq K\\ (forced) | \\(y - \mathrm{lim}) / K\\  | \\1 - \epsilon\\            |
| 3          | type `"m"`            | \\(y - \mathrm{lim}) / K\\  | \\(y + \mathrm{lim}) / K\\  |
| 3          | type `"l"`            | \\(y - 2\mathrm{lim}) / K\\ | \\y / K\\                   |
| 3          | type `"r"`            | \\y / K\\                   | \\(y + 2\mathrm{lim}) / K\\ |

All endpoints are clamped to \\\[\epsilon, 1 - \epsilon\]\\ with
\\\epsilon = 10^{-5}\\ to avoid boundary issues in the beta likelihood.

The midpoint approximation `yt` is computed as:

- \\y_t = y\\ when \\y \in (0, 1)\\ (continuous data).

- \\y_t = y / K\\ when \\y\\ is on the integer scale.

This value is used exclusively as an initialization aid for
starting-value computation and does not enter the likelihood.

**Interaction with the fitting pipeline**:

This function is called internally by `.extract_response()` when the
data does *not* carry the `"bs_prepared"` attribute. When data has been
pre-processed by
[`bs_prepare`](https://evandeilton.github.io/betaregscale/reference/bs_prepare.md)
or by simulation with forced delta
([`betaregscale_simulate`](https://evandeilton.github.io/betaregscale/reference/betaregscale_simulate.md)
with `delta != NULL`), the pre-computed columns are used directly and
`check_response()` is skipped.

## See also

[`bs_prepare`](https://evandeilton.github.io/betaregscale/reference/bs_prepare.md)
for the analyst-facing pre-processing function;
[`betaregscale_simulate`](https://evandeilton.github.io/betaregscale/reference/betaregscale_simulate.md)
for simulation with forced delta.

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

# Force delta = 1 on non-boundary observations:
# endpoints use actual y values, preserving variation
y2 <- c(30, 60)
check_response(y2, ncuts = 100, delta = c(1L, 1L))
#>       left right  yt  y delta
#> [1,] 1e-05 0.305 0.3 30     1
#> [2,] 1e-05 0.605 0.6 60     1
#  left = (eps, eps), right = (30.5/100, 60.5/100)
```
