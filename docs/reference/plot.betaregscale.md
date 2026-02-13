# Diagnostic plots for beta interval regression

Produces up to six diagnostic plots for a fitted `"betaregscale"` model:
residuals vs indices, Cook's distance, residuals vs linear predictor,
residuals vs fitted values, a half-normal plot with simulated envelope,
and predicted vs observed.

## Usage

``` r
# S3 method for class 'betaregscale'
plot(
  x,
  which = 1:4,
  type = "rqr",
  nsim = 100L,
  level = 0.9,
  caption = c("Residuals vs indices", "Cook's distance", "Residuals vs linear predictor",
    "Residuals vs fitted values", "Half-normal plot", "Predicted vs observed"),
  sub.caption = NULL,
  ask = prod(par("mfcol")) < length(which) && dev.interactive(),
  gg = FALSE,
  ...
)
```

## Arguments

- x:

  A fitted `"betaregscale"` object.

- which:

  Integer vector selecting which plots to draw (default `1:4`).

- type:

  Character: residual type passed to
  [`residuals.betaregscale`](https://evandeilton.github.io/betaregscale/reference/residuals.betaregscale.md)
  (default `"rqr"`).

- nsim:

  Integer: number of simulations for the half-normal envelope (default
  100).

- level:

  Numeric: confidence level for the envelope (default 0.9).

- caption:

  Character vector of panel captions.

- sub.caption:

  Subtitle; defaults to the model call.

- ask:

  Logical: prompt before each page of plots?

- gg:

  Logical: use ggplot2? (default `FALSE`).

- ...:

  Further arguments passed to base
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly returns `x`.
