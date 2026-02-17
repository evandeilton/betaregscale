# Diagnostic plots for mixed beta interval regression

Produces diagnostic plots for fitted `"brsmm"` models: residuals vs
indices, Cook's distance, residuals vs linear predictor, residuals vs
fitted values, half-normal envelope, and predicted vs observed.

## Usage

``` r
# S3 method for class 'brsmm'
plot(
  x,
  which = 1:4,
  type = c("response", "pearson"),
  nsim = 100L,
  level = 0.9,
  caption = c("Residuals vs indices", "Cook's distance", "Residuals vs linear predictor",
    "Residuals vs fitted values", "Half-normal plot", "Predicted vs observed"),
  sub.caption = NULL,
  ask = prod(par("mfcol")) < length(which) && dev.interactive(),
  gg = FALSE,
  title = NULL,
  theme = NULL,
  ...
)
```

## Arguments

- x:

  A fitted `"brsmm"` object.

- which:

  Integer vector selecting which panels to draw (default `1:4`).

- type:

  Residual type passed to
  [`residuals.brsmm`](https://evandeilton.github.io/betaregscale/reference/residuals.brsmm.md)
  (`"response"` or `"pearson"`).

- nsim:

  Number of simulations for half-normal envelope.

- level:

  Confidence level for the half-normal envelope.

- caption:

  Character vector of plot captions.

- sub.caption:

  Optional subtitle; defaults to model call.

- ask:

  Logical: prompt before each new page?

- gg:

  Logical: use ggplot2 backend?

- title:

  Optional global title for ggplot output. If `NULL`, panel captions are
  used.

- theme:

  Optional ggplot2 theme object (e.g.,
  [`ggplot2::theme_bw()`](https://ggplot2.tidyverse.org/reference/ggtheme.html)).
  If `NULL`, a minimal theme is used.

- ...:

  Further arguments passed to base
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly returns `x`.
