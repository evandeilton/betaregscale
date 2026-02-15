# Graphical and tabular censoring summary

Produces a visual summary of the censoring structure in a fitted `"brs"`
model or a response matrix produced by
[`brs_check`](https://evandeilton.github.io/betaregscale/reference/brs_check.md).
The summary includes:

1.  Bar chart of censoring type counts

2.  Histogram of midpoint responses colored by censoring type

3.  Interval plot showing \\\[l_i, u_i\]\\ segments

4.  Proportion table of censoring types

## Usage

``` r
brs_cens(object, n_sample = 100L, gg = FALSE, ...)
```

## Arguments

- object:

  A fitted `"brs"` object, a matrix returned by
  [`brs_check`](https://evandeilton.github.io/betaregscale/reference/brs_check.md),
  or a data frame returned by
  [`brs_prep`](https://evandeilton.github.io/betaregscale/reference/brs_prep.md)
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
y <- c(0, 3, 5, 7, 10)
Y <- brs_check(y, ncuts = 10)
brs_cens(Y)


prep <- brs_prep(data.frame(y = y), ncuts = 10)
#> brs_prep: n = 5 | exact = 0, left = 1, right = 1, interval = 3
brs_cens(prep)
```
