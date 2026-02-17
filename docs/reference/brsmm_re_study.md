# Random-effects study for brsmm models

Provides a compact numeric study of random effects, including: estimated
covariance matrix, correlation matrix, per-term standard deviations,
empirical mean/SD of posterior modes, shrinkage ratio, and a normality
check by Shapiro-Wilk (when applicable).

## Usage

``` r
brsmm_re_study(object, ...)
```

## Arguments

- object:

  A fitted `"brsmm"` object.

- ...:

  Currently ignored.

## Value

A list with class `"brsmm_re_study"`.

## Examples

``` r
# \donttest{
set.seed(123)
g <- 10
ni <- 8
id <- factor(rep(seq_len(g), each = ni))
n <- length(id)
x1 <- rnorm(n)
b0 <- rnorm(g, sd = 0.4)
b1 <- rnorm(g, sd = 0.2)
mu <- plogis(0.2 + 0.5 * x1 + b0[id] + b1[id] * x1)
phi <- plogis(-0.2)
shp <- brs_repar(mu = mu, phi = rep(phi, n), repar = 2)
y <- round(stats::rbeta(n, shp$shape1, shp$shape2) * 100)
d <- data.frame(y = y, x1 = x1, id = id)
fit <- brsmm(y ~ x1, random = ~ 1 + x1 | id, data = d, repar = 2)

rs <- brsmm_re_study(fit)
print(rs)
#> 
#> Random-effects study
#> Groups: 10 
#> 
#> Summary by term:
#>         term sd_model mean_mode sd_mode shrinkage_ratio shapiro_p
#>  (Intercept)   0.0014         0       0               0    0.9957
#>           x1   0.0032         0       0               0    0.7480
#> 
#> Estimated covariance matrix D:
#>      [,1] [,2]
#> [1,]    0    0
#> [2,]    0    0
#> 
#> Estimated correlation matrix:
#>        [,1]   [,2]
#> [1,] 1.0000 0.2508
#> [2,] 0.2508 1.0000
knitr::kable(rs$summary, digits = 4)
#> 
#> 
#> |term        | sd_model| mean_mode| sd_mode| shrinkage_ratio| shapiro_p|
#> |:-----------|--------:|---------:|-------:|---------------:|---------:|
#> |(Intercept) |   0.0014|         0|       0|               0|    0.9957|
#> |x1          |   0.0032|         0|       0|               0|    0.7480|
# }
```
