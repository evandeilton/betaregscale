# Coefficient estimates with inference

Coefficient estimates with inference

## Usage

``` r
brs_est(object, alpha = 0.05)
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
sim <- brs_sim(
  formula = ~x1, data = data.frame(x1 = rnorm(50)),
  beta = c(0, 0.5), phi = 0.1, ncuts = 10, repar = 2
)
fit <- brs(y ~ x1, data = sim, repar = 2)
brs_est(fit)
#>      variable   estimate        se   z_value      p_value    ci_lower
#> 1 (Intercept) -2.8611765 0.1361136 -21.02050 4.258938e-98 -3.12795428
#> 2          x1  0.2317141 0.1088125   2.12948 3.321453e-02  0.01844551
#> 3       (phi) -3.0117569 0.2372581 -12.69401 6.383324e-37 -3.47677434
#>     ci_upper
#> 1 -2.5943987
#> 2  0.4449827
#> 3 -2.5467395
# }
```
