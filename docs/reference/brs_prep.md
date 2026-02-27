# Pre-process analyst data for beta interval regression

Validates and transforms raw data into the format required by
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md).
The analyst can supply data in several ways:

1.  **Minimal (Mode 1)**: only the score `y`. Censoring is inferred
    automatically: \\y = 0 \to \delta = 1\\, \\y = K \to \delta = 2\\,
    \\0 \< y \< K \to \delta = 3\\, \\y \in (0, 1) \to \delta = 0\\.

2.  **Classic (Mode 2)**: `y` + explicit `delta`. The analyst declares
    the censoring type; interval endpoints are computed using the actual
    `y` value.

3.  **Interval (Mode 3)**: `left` and/or `right` columns (on the
    original scale). Censoring is inferred from the NA pattern.

4.  **Full (Mode 4)**: `y`, `left`, and `right` together. The analyst's
    own endpoints are rescaled directly to \\(0, 1)\\.

All covariate columns are preserved unchanged in the output.

## Usage

``` r
brs_prep(
  data,
  y = "y",
  delta = "delta",
  left = "left",
  right = "right",
  ncuts = 100L,
  lim = 0.5
)
```

## Arguments

- data:

  Data frame containing the response variable and covariates.

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
`"is_prepared"` (`TRUE`), `"ncuts"` and `"lim"` so that
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md) can
detect prepared data and skip the internal
[`brs_check`](https://evandeilton.github.io/betaregscale/reference/brs_check.md)
call.

## Details

**Priority rule**: if `delta` is provided (non-`NA`), it takes
precedence over all automatic classification rules. When `delta` is
`NA`, the function infers the censoring type from the pattern of `left`,
`right`, and `y`:

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

**Endpoint formulas for Mode 2 (y + explicit delta)**:

When the analyst supplies `delta` explicitly, the endpoint computation
uses the actual `y` value to produce observation-specific bounds. This
is the same logic used by
[`brs_check`](https://evandeilton.github.io/betaregscale/reference/brs_check.md)
with a user-supplied `delta` vector:

|  |  |  |  |
|----|----|----|----|
| \\\delta\\ | Condition | \\l_i\\ (left) | \\u_i\\ (right) |
| 0 | (any) | \\y / K\\ | \\y / K\\ |
| 1 | \\y = 0\\ | \\\epsilon\\ | \\\mathrm{lim} / K\\ |
| 1 | \\y \neq 0\\ | \\\epsilon\\ | \\(y + \mathrm{lim}) / K\\ |
| 2 | \\y = K\\ | \\(K - \mathrm{lim}) / K\\ | \\1 - \epsilon\\ |
| 2 | \\y \neq K\\ | \\(y - \mathrm{lim}) / K\\ | \\1 - \epsilon\\ |
| 3 | type `"m"` | \\(y - \mathrm{lim}) / K\\ | \\(y + \mathrm{lim}) / K\\ |

**Consistency warnings**: when the analyst supplies `delta` values that
are unusual for the given `y` (e.g., \\\delta = 1\\ but \\y \neq 0\\),
the function emits a warning but proceeds. This is by design for Monte
Carlo workflows where forced delta on non-boundary observations is
intentional.

All endpoints are clamped to \\\[\epsilon, 1 - \epsilon\]\\ with
\\\epsilon = 10^{-5}\\.

## References

Lopes, J. E. (2023). *Modelos de regressao beta para dados de escala*.
Master's dissertation, Universidade Federal do Parana, Curitiba. URI:
<https://hdl.handle.net/1884/86624>.

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

[`brs_check`](https://evandeilton.github.io/betaregscale/reference/brs_check.md)
for the automatic classification of raw scale scores;
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md) for
fitting the model.

## Examples

``` r
# --- Mode 1: y only (automatic classification, like brs_check) ---
d1 <- data.frame(y = c(0, 3, 5, 7, 10), x1 = rnorm(5))
brs_prep(d1, ncuts = 10)
#> brs_prep: n = 5 | exact = 0, left = 1, right = 1, interval = 3
#>      left   right      yt  y delta          x1
#> 1 0.00001 0.05000 0.00001  0     1 -1.25150957
#> 2 0.25000 0.35000 0.30000  3     3  0.52848796
#> 3 0.45000 0.55000 0.50000  5     3 -1.24761627
#> 4 0.65000 0.75000 0.70000  7     3 -0.04165134
#> 5 0.95000 0.99999 0.99999 10     2 -1.05473729

# --- Mode 2: y + explicit delta ---
d2 <- data.frame(
  y = d1$y,
  delta = c(0, 3, 3, 3, 0), # Force interval-censoring for 3,5,7
  x1 = d1$x1
)
brs_prep(d2, ncuts = 100)
#> brs_prep: n = 5 | exact = 2, left = 0, right = 0, interval = 3
#>      left   right    yt  y delta          x1
#> 1 0.00001 0.00001 1e-05  0     0 -1.25150957
#> 2 0.02500 0.03500 3e-02  3     3  0.52848796
#> 3 0.04500 0.05500 5e-02  5     3 -1.24761627
#> 4 0.06500 0.07500 7e-02  7     3 -0.04165134
#> 5 0.10000 0.10000 1e-01 10     0 -1.05473729

# --- Mode 3: left/right with NA patterns ---
d3 <- data.frame(
  left = c(NA, 20, 30, NA),
  right = c(5, NA, 45, NA),
  y = c(NA, NA, NA, 50),
  x1 = d1$x1[1:4]
)
brs_prep(d3, ncuts = 100)
#> brs_prep: n = 4 | exact = 1, left = 1, right = 1, interval = 1
#>    left   right    yt  y delta          x1
#> 1 1e-05 0.05000 0.025 NA     1 -1.25150957
#> 2 2e-01 0.99999 0.600 NA     2  0.52848796
#> 3 3e-01 0.45000 0.375 NA     3 -1.24761627
#> 4 5e-01 0.50000 0.500 50     0 -0.04165134

# --- Mode 4: y + left + right (analyst-supplied intervals) ---
d4 <- data.frame(
  y = c(50, 75),
  left = c(48, 73),
  right = c(52, 77),
  x1 = rnorm(2)
)
brs_prep(d4, ncuts = 100)
#> brs_prep: n = 2 | exact = 0, left = 0, right = 0, interval = 2
#>   left right   yt  y delta          x1
#> 1 0.48  0.52 0.50 50     3 -0.09817368
#> 2 0.73  0.77 0.75 75     3  0.26143179

# --- Fitting after prep ---
# \donttest{
dat5 <- data.frame(
  y = c(
    0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
    10, 40, 55, 70, 85, 25, 35, 65, 80, 15
  ),
  x1 = rep(c(1, 2), 10)
)
prep5 <- brs_prep(dat5, ncuts = 100)
#> brs_prep: n = 20 | exact = 0, left = 1, right = 1, interval = 18
fit5 <- brs(y ~ x1, data = prep5)
summary(fit5)
#> 
#> Call:
#> brs(formula = y ~ x1, data = prep5)
#> 
#> Quantile residuals:
#>     Min      1Q  Median      3Q     Max 
#> -2.2706 -0.4813  0.0555  0.5455  3.1621 
#> 
#> Coefficients (mean model with logit link):
#>             Estimate Std. Error z value Pr(>|z|)
#> (Intercept)   0.2551     0.8644   0.295    0.768
#> x1           -0.2202     0.5412  -0.407    0.684
#> 
#> Phi coefficients (precision model with logit link):
#>       Estimate Std. Error z value Pr(>|z|)
#> (phi)  -0.3929     0.2763  -1.422    0.155
#> ---
#> Log-likelihood: -92.6521 on 3 Df | AIC: 191.3041 | BIC: 194.2913 
#> Pseudo R-squared: 0.0029 
#> Number of iterations: 17 (BFGS) 
#> Censoring: 18 interval | 1 left | 1 right 
#> 
# }
```
