# Graphical and tabular censoring summary

Produces a visual summary of the censoring structure in a fitted
`"betaregscale"` model or a response matrix produced by
[`check_response`](https://evandeilton.github.io/betaregscale/reference/check_response.md).
The summary includes:

1.  Bar chart of censoring type counts

2.  Histogram of midpoint responses colored by censoring type

3.  Interval plot showing \\\[l_i, u_i\]\\ segments

4.  Proportion table of censoring types

## Usage

``` r
censoring_summary(object, n_sample = 100L, gg = FALSE, ...)
```

## Arguments

- object:

  A fitted `"betaregscale"` object, a matrix returned by
  [`check_response`](https://evandeilton.github.io/betaregscale/reference/check_response.md),
  or a data frame returned by
  [`bs_prepare`](https://evandeilton.github.io/betaregscale/reference/bs_prepare.md)
  (must contain columns `left`, `right`, `yt`, and `delta`).

- n_sample:

  Integer: maximum number of observations to show in the interval plot
  (default 100). If the data has more observations, a random sample is
  drawn.

- gg:

  Logical: use ggplot2? (default `FALSE`).

- ...:

  Further arguments (currently ignored).

## Value

Invisibly returns a data frame with censoring counts and proportions.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- betaregscale(y ~ x1 + x2, data = sim)
censoring_summary(fit)
censoring_summary(fit, gg = TRUE)
} # }
```
