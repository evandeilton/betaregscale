# Pre-process analyst data for beta interval regression

Validates and transforms raw data into the format required by
[`betaregscale`](https://evandeilton.github.io/betaregscale/reference/betaregscale.md).
The analyst can supply data in several ways:

1.  **Minimal**: only the score `y`. Censoring is inferred automatically
    (equivalent to
    [`check_response`](https://evandeilton.github.io/betaregscale/reference/check_response.md)).

2.  **Classic**: `y` + explicit `delta`. The analyst states the
    censoring type; interval endpoints are computed.

3.  **Interval**: `left` and/or `right` columns (on the original scale).
    Censoring is inferred from the NA pattern.

4.  **Full**: `y`, `left`, and `right` together. The analyst's own
    endpoints are rescaled directly to \\(0, 1)\\.

All covariate columns are preserved unchanged in the output data frame.

## Usage

``` r
bs_prepare(
  data,
  y = "y",
  delta = "delta",
  left = "left",
  right = "right",
  ncuts = 100L,
  type = "m",
  lim = 0.5
)
```

## Arguments

- data:

  A `data.frame` containing the response and (optionally) covariates.

- y:

  Character: name of the score column (default `"y"`).

- delta:

  Character: name of the censoring indicator column (default `"delta"`).
  Values must be in `{0, 1, 2, 3}`.

- left:

  Character: name of the left-endpoint column (default `"left"`).

- right:

  Character: name of the right-endpoint column (default `"right"`).

- ncuts:

  Integer: number of scale categories (default 100).

- type:

  Character: interval type for interior scores when only `y` and `delta`
  are available. `"m"` = midpoint (default), `"l"` = left-aligned, `"r"`
  = right-aligned. See
  [`check_response`](https://evandeilton.github.io/betaregscale/reference/check_response.md).

- lim:

  Numeric: half-width of the uncertainty region (default 0.5). Used only
  when constructing intervals from `y` alone.

## Value

A `data.frame` with the following columns appended or replaced:

- `left`:

  Lower endpoint on \\(0, 1)\\.

- `right`:

  Upper endpoint on \\(0, 1)\\.

- `yt`:

  Midpoint approximation on \\(0, 1)\\.

- `y`:

  Original scale value (preserved for reference).

- `delta`:

  Censoring indicator: 0 = exact, 1 = left, 2 = right, 3 = interval.

Covariate columns are preserved. The output carries attributes
`"bs_prepared"` (`TRUE`), `"ncuts"`, `"type"`, and `"lim"` so that
[`betaregscale`](https://evandeilton.github.io/betaregscale/reference/betaregscale.md)
can detect prepared data and skip the internal
[`check_response`](https://evandeilton.github.io/betaregscale/reference/check_response.md)
call.

## Details

**Priority rule**: if `delta` is provided (non-`NA`), it takes
precedence. When `delta` is `NA`, the function infers the censoring type
from the pattern of `left`, `right`, and `y`:

|        |         |      |         |                              |                     |
|--------|---------|------|---------|------------------------------|---------------------|
| `left` | `right` | `y`  | `delta` | Interpretation               | Inferred \\\delta\\ |
| `NA`   | 5       | `NA` | `NA`    | Left-censored (below 5)      | 1                   |
| 20     | `NA`    | `NA` | `NA`    | Right-censored (above 20)    | 2                   |
| 30     | 45      | `NA` | `NA`    | Interval-censored \[30, 45\] | 3                   |
| `NA`   | `NA`    | 50   | `NA`    | Exact observation            | 0                   |
| `NA`   | `NA`    | 50   | 3       | Analyst says interval        | 3                   |
| `NA`   | `NA`    | 0    | 1       | Analyst says left-censored   | 1                   |
| `NA`   | `NA`    | 99   | 2       | Analyst says right-censored  | 2                   |

When `y`, `left`, and `right` are all present for the same observation,
the analyst's `left`/`right` values are used directly (rescaled by \\K
=\\ `ncuts`) and `delta` is set to 3 (interval-censored) unless the
analyst supplied `delta` explicitly.

All endpoints are clamped to \\\[\epsilon, 1 - \epsilon\]\\ with
\\\epsilon = 10^{-5}\\.

## See also

[`check_response`](https://evandeilton.github.io/betaregscale/reference/check_response.md)
for the automatic classification of raw scale scores;
[`betaregscale`](https://evandeilton.github.io/betaregscale/reference/betaregscale.md)
for fitting the model.

## Examples

``` r
# --- Mode 1: y only (automatic classification, like check_response) ---
d1 <- data.frame(y = c(0, 3, 5, 7, 10), x1 = rnorm(5))
bs_prepare(d1, ncuts = 10)
#> bs_prepare: n = 5 | exact = 0, left = 1, right = 1, interval = 3
#>      left   right      yt  y delta         x1
#> 1 0.00001 0.05000 0.00001  0     1  0.9688992
#> 2 0.25000 0.35000 0.30000  3     3 -0.8757371
#> 3 0.45000 0.55000 0.50000  5     3 -2.1360245
#> 4 0.65000 0.75000 0.70000  7     3 -1.5228276
#> 5 0.95000 0.99999 0.99999 10     2 -1.1131182

# --- Mode 2: y + explicit delta ---
d2 <- data.frame(
  y     = c(50, 0, 99, 50),
  delta = c( 0, 1,  2,  3),
  x1    = rnorm(4)
)
bs_prepare(d2, ncuts = 100)
#> Warning: Observation(s) 3: delta = 2 (right-censored) but y != 100.
#> bs_prepare: n = 4 | exact = 1, left = 1, right = 1, interval = 1
#>      left   right      yt  y delta           x1
#> 1 0.50000 0.50000 0.50000 50     0  1.240983896
#> 2 0.00001 0.00500 0.00001  0     1  0.003481935
#> 3 0.99500 0.99999 0.99999 99     2 -1.237794586
#> 4 0.49500 0.50500 0.50000 50     3  0.555699119

# --- Mode 3: left/right with NA patterns ---
d3 <- data.frame(
  left  = c(NA, 20, 30, NA),
  right = c( 5, NA, 45, NA),
  y     = c(NA, NA, NA, 50),
  x1    = rnorm(4)
)
bs_prepare(d3, ncuts = 100)
#> bs_prepare: n = 4 | exact = 0, left = 1, right = 1, interval = 2
#>      left   right    yt  y delta         x1
#> 1 0.00001 0.05000 0.025 NA     1 -2.1831493
#> 2 0.20000 0.99999 0.600 NA     2 -0.2470245
#> 3 0.30000 0.45000 0.375 NA     3  1.1128569
#> 4 0.49500 0.50500 0.500 50     3 -0.3416732

# --- Mode 4: y + left + right (analyst-supplied intervals) ---
d4 <- data.frame(
  y     = c(50, 75),
  left  = c(48, 73),
  right = c(52, 77),
  x1    = rnorm(2)
)
bs_prepare(d4, ncuts = 100)
#> bs_prepare: n = 2 | exact = 0, left = 0, right = 0, interval = 2
#>   left right   yt  y delta       x1
#> 1 0.48  0.52 0.50 50     3 1.305438
#> 2 0.73  0.77 0.75 75     3 1.795324

# --- End-to-end workflow ---
# \donttest{
set.seed(42)
n <- 200
dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
sim <- betaregscale_simulate(
  formula = ~ x1 + x2, data = dat,
  beta = c(0.2, -0.5, 0.3), phi = 1 / 5
)
prep <- bs_prepare(sim, ncuts = 100)
#> bs_prepare: n = 200 | exact = 0, left = 18, right = 21, interval = 161
fit <- betaregscale(y ~ x1 + x2, data = prep)
summary(fit)
#> 
#> Call:
#> betaregscale(formula = y ~ x1 + x2, data = prep)
#> 
#> Quantile residuals:
#>     Min      1Q  Median      3Q     Max 
#> -3.0625 -0.5896  0.2555  0.6723  1.5528 
#> 
#> Coefficients (mean model with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -5.15677    0.08757 -58.886  < 2e-16 ***
#> x1          -0.38266    0.06216  -6.156 7.45e-10 ***
#> x2           0.12275    0.06556   1.872   0.0612 .  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Phi coefficients (precision model with logit link):
#>       Estimate Std. Error z value Pr(>|z|)    
#> (phi)  -4.8810     0.1352  -36.11   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> ---
#> Log-likelihood: -896.1340 on 4 Df
#> Pseudo R-squared: 0.1097 
#> Number of iterations: 29 (BFGS) 
#> Censoring: 161 interval | 18 left | 21 right 
#> 
# }
```
