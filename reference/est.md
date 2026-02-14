# Coefficient estimates with inference

Coefficient estimates with inference

## Usage

``` r
est(object, alpha = 0.05)
```

## Arguments

- object:

  A fitted `"betaregscale"` object.

- alpha:

  Significance level (default 0.05).

## Value

Data frame of estimates, standard errors, z-values, and p-values.

## Examples

``` r
# \donttest{
sim <- betaregscale_simulate(
  formula = ~x1, data = data.frame(x1 = rnorm(50)),
  beta = c(0, 0.5), phi = 0.1, ncuts = 10, repar = 2
)
fit <- betaregscale(y ~ x1, data = sim, repar = 2)
#> Warning: The 'type' argument of betaregscale_fit() is deprecated and will be removed in a future version. Use bs_prepare() to control interval geometry.
est(fit)
#>      variable  estimate        se    z_value       p_value    ci_lower
#> 1 (Intercept) -2.923978 0.1321935 -22.118922 2.078346e-108 -3.18307302
#> 2          x1  0.334766 0.1309703   2.556045  1.058694e-02  0.07806892
#> 3       (phi) -3.150176 0.2370857 -13.287078  2.751150e-40 -3.61485552
#>     ci_upper
#> 1 -2.6648839
#> 2  0.5914631
#> 3 -2.6854967
# }
```
