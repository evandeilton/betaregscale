# Extract the Hessian matrix

Extract the Hessian matrix

## Usage

``` r
brs_hessian(object)
```

## Arguments

- object:

  A fitted `"betaregscale"` object.

## Value

Numeric Hessian matrix.

## Examples

``` r
# \donttest{
sim <- brs_sim(
  formula = ~x1, data = data.frame(x1 = rnorm(50)),
  beta = c(0, 0.5), phi = 0.1, ncuts = 10, repar = 2
)
fit <- brs(y ~ x1, data = sim, repar = 2)
brs_hessian(fit)
#>             (Intercept)          x1      (phi)
#> (Intercept) -67.1756949  -0.7715532  20.695689
#> x1           -0.7715532 -72.3879779  -2.322661
#> (phi)        20.6956894  -2.3226607 -23.549628
# }
```
