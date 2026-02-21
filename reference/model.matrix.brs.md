# Extract design matrix

Extract design matrix

## Usage

``` r
# S3 method for class 'brs'
model.matrix(object, model = c("mean", "precision"), ...)
```

## Arguments

- object:

  A fitted `"betaregscale"` object.

- model:

  Character: `"mean"` (default) or `"precision"`.

- ...:

  Ignored.

## Value

The design matrix for the specified submodel.

## See also

[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md),
[`formula.brs`](https://evandeilton.github.io/betaregscale/reference/formula.brs.md),
[`coef.brs`](https://evandeilton.github.io/betaregscale/reference/coef.brs.md)

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
head(model.matrix(fit))
#>   (Intercept) x1
#> 1           1  1
#> 2           1  2
#> 3           1  1
#> 4           1  2
#> 5           1  1
#> 6           1  2
head(model.matrix(fit, model = "precision"))
#>      (Intercept)
#> [1,]           1
#> [2,]           1
#> [3,]           1
#> [4,]           1
#> [5,]           1
#> [6,]           1
# }
```
