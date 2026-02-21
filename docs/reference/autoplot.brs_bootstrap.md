# ggplot2 autoplot for bootstrap results

Produces visual summaries for objects returned by
[`brs_bootstrap`](https://evandeilton.github.io/betaregscale/reference/brs_bootstrap.md).

## Usage

``` r
# S3 method for class 'brs_bootstrap'
autoplot(
  object,
  type = c("ci_forest", "dist", "qq", "stability"),
  parameter = NULL,
  title = NULL,
  caption = NULL,
  max_parameters = 12L,
  ci_level = NULL,
  theme = NULL,
  ...
)
```

## Arguments

- object:

  An object of class `"brs_bootstrap"`.

- type:

  Plot type: `"ci_forest"`, `"dist"`, `"qq"`, or `"stability"`.

- parameter:

  Optional parameter name used by `type = "dist"`, `"qq"`, and
  `"stability"`. If `NULL`, the first parameter is used.

- title:

  Optional plot title override.

- caption:

  Optional subtitles/titles for plot types. Accepts:

  - a single string (used for the selected `type`);

  - a character vector/list with up to four entries in the order
    `ci_forest`, `dist`, `qq`, `stability`.

- max_parameters:

  Maximum number of parameters shown in `type = "ci_forest"`.

- ci_level:

  Confidence level used in `type = "stability"`. Defaults to the level
  stored in `object`.

- theme:

  Optional ggplot2 theme object (e.g.,
  [`ggplot2::theme_bw()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)).
  If `NULL`,
  [`ggplot2::theme_minimal()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)
  is used.

- ...:

  Currently ignored.

## Value

A `ggplot2` object.

## Details

For `type = "dist"`, `"qq"`, and `"stability"`, bootstrap draws must be
present in `attr(object, "boot_draws")`, obtained by fitting with
`brs_bootstrap(..., keep_draws = TRUE)`.

## See also

[`brs_bootstrap`](https://evandeilton.github.io/betaregscale/reference/brs_bootstrap.md),
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md),
[`autoplot.brs`](https://evandeilton.github.io/betaregscale/reference/autoplot.brs.md)

## Examples

``` r
# \donttest{
dat <- data.frame(
  y = c(
    0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
    10, 40, 55, 70, 85, 25, 35, 65, 80, 15
  ),
  x1 = rep(c(1, 2), 10)
)
prep <- brs_prep(dat, ncuts = 100)
#> brs_prep: n = 20 | exact = 0, left = 1, right = 1, interval = 18
fit <- brs(y ~ x1, data = prep)
boot <- brs_bootstrap(fit, R = 50)
ggplot2::autoplot(boot, type = "ci_forest")

# }
```
