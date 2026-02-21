# Extract design matrix

Extract design matrix

## Usage

``` r
# S3 method for class 'brsmm'
model.matrix(object, model = c("mean", "precision", "random"), ...)
```

## Arguments

- object:

  A fitted `"brsmm"` object.

- model:

  Character: `"mean"` (default), `"precision"`, or `"random"`.

- ...:

  Ignored.

## Value

The design matrix for the specified submodel.

## See also

[`brsmm`](https://evandeilton.github.io/betaregscale/reference/brsmm.md),
[`formula.brsmm`](https://evandeilton.github.io/betaregscale/reference/formula.brsmm.md)

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
fit <- brsmm(y ~ x1, random = ~ 1 | id, data = prep)
head(model.matrix(fit))
#>   (Intercept) x1
#> 1           1  1
#> 2           1  2
#> 3           1  1
#> 4           1  2
#> 5           1  1
#> 6           1  2
head(model.matrix(fit, model = "random"))
#>   (Intercept)
#> 1           1
#> 2           1
#> 3           1
#> 4           1
#> 5           1
#> 6           1
# }
```
