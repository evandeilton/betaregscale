# ggplot2 autoplot for marginal effects

Produces visual summaries for objects returned by
[`brs_marginaleffects`](https://evandeilton.github.io/betaregscale/reference/brs_marginaleffects.md).

## Usage

``` r
# S3 method for class 'brs_marginaleffects'
autoplot(
  object,
  type = c("forest", "magnitude", "dist"),
  variable = NULL,
  top_n = 12L,
  title = NULL,
  caption = NULL,
  theme = NULL,
  ...
)
```

## Arguments

- object:

  An object of class `"brs_marginaleffects"`.

- type:

  Plot type: `"forest"`, `"magnitude"`, or `"dist"`.

- variable:

  Optional variable name for `type = "dist"`.

- top_n:

  Maximum number of variables shown in `"magnitude"` (ordered by
  `|AME|`).

- title:

  Optional plot title override.

- caption:

  Optional subtitle override.

- theme:

  Optional ggplot2 theme object. If `NULL`,
  [`ggplot2::theme_minimal()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)
  is used.

- ...:

  Currently ignored.

## Value

A `ggplot2` object.

## Details

`type = "dist"` requires AME simulation draws stored in
`attr(object, "ame_draws")`, which are available when marginal effects
are computed with `keep_draws = TRUE` and `interval = TRUE`.
