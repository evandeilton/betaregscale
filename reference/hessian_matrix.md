# Extract the Hessian matrix

Extract the Hessian matrix

## Usage

``` r
hessian_matrix(object)
```

## Arguments

- object:

  A fitted `"betaregscale"` object.

## Value

Numeric Hessian matrix.

## Examples

``` r
# \donttest{
sim <- betaregscale_simulate(
  formula = ~x1, data = data.frame(x1 = rnorm(50)),
  beta = c(0, 0.5), phi = 0.1, ncuts = 10, repar = 2
)
fit <- betaregscale(y ~ x1, data = sim, repar = 2)
#> Warning: The 'type' argument of betaregscale_fit() is deprecated and will be removed in a future version. Use bs_prepare() to control interval geometry.
hessian_matrix(fit)
#>             (Intercept)         x1      (phi)
#> (Intercept)   -63.05945 -35.787161  20.709762
#> x1            -35.78716 -97.998177   8.009939
#> (phi)          20.70976   8.009939 -23.473559
# }
```
