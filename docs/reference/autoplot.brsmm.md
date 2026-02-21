# ggplot2 autoplot for brsmm models

Produces ggplot2 diagnostics tailored to mixed beta interval models.

## Usage

``` r
# S3 method for class 'brsmm'
autoplot(
  object,
  type = c("calibration", "score_dist", "ranef_qq", "residuals_by_group",
    "ranef_caterpillar", "ranef_density", "ranef_pairs"),
  bins = 10L,
  scores = NULL,
  residual_type = c("response", "pearson"),
  max_groups = 25L,
  ...
)
```

## Arguments

- object:

  A fitted `"brsmm"` object.

- type:

  Plot type: `"calibration"`, `"score_dist"`, `"ranef_qq"`, or
  `"residuals_by_group"`, `"ranef_caterpillar"`, `"ranef_density"`,
  `"ranef_pairs"`.

- bins:

  Number of bins used in calibration plots.

- scores:

  Optional integer vector of scores for `"score_dist"`. Defaults to all
  scores from `0` to `ncuts`.

- residual_type:

  Residual type passed to
  [`residuals.brsmm`](https://evandeilton.github.io/betaregscale/reference/residuals.brsmm.md)
  for `type = "residuals_by_group"`.

- max_groups:

  Maximum number of groups displayed in `"residuals_by_group"`.

- ...:

  Currently ignored.

## Value

A `ggplot2` object.

## Examples

``` r
# \donttest{
dat <- data.frame(
  y = c(
    0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
    10, 40, 55, 70, 85, 25, 35, 65, 80, 15
  ),
  x1 = rep(c(1, 2), 10),
  id = factor(rep(1:4, each = 5))
)
prep <- brs_prep(dat, ncuts = 100)
#> brs_prep: n = 20 | exact = 0, left = 1, right = 1, interval = 18
fit_mm <- brsmm(y ~ x1, random = ~ 1 | id, data = prep)
ggplot2::autoplot(fit_mm, type = "calibration", bins = 4)

ggplot2::autoplot(fit_mm, type = "score_dist")

ggplot2::autoplot(fit_mm, type = "ranef_qq")

ggplot2::autoplot(fit_mm, type = "ranef_caterpillar")

ggplot2::autoplot(fit_mm, type = "ranef_density")

# }
```
