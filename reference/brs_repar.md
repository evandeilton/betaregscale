# Beta shape parameters from the parameters of each scheme

Converts the pair `(mu, phi)` of one of three parameterisations into the
shapes \\(a, b)\\ of the beta density \\f(y) = y^{a - 1}(1 - y)^{b -
1}/B(a, b)\\, \\0 \< y \< 1\\.

## Usage

``` r
brs_repar(mu, phi, repar = 2L)
```

## Arguments

- mu:

  Numeric vector: the first parameter, i.e. the mean in \\(0, 1)\\ for
  `repar = 1, 2` and the shape \\p \> 0\\ for `repar = 0`.

- phi:

  Numeric vector (or scalar): the second parameter, i.e. the precision
  \\\phi \> 0\\, the dispersion \\\phi \in (0, 1)\\ or the shape \\q \>
  0\\.

- repar:

  Integer (0, 1 or 2) selecting the scheme (default 2).

## Value

A `data.frame` with columns `shape1` (\\a\\) and `shape2` (\\b\\), one
row per element of the recycled inputs.

## Details

|  |  |  |  |  |
|----|----|----|----|----|
| `repar` | `mu`, `phi` | \\(a, b)\\ | \\E\[Y\]\\ | \\\mathrm{Var}\[Y\]\\ |
| 0 | shapes \\p, q \> 0\\ | \\(p, q)\\ | \\p/(p + q)\\ | \\pq/\\(p + q)^2 (p + q + 1)\\\\ |
| 1 | mean \\\mu \in (0, 1)\\, precision \\\phi \> 0\\ | \\(\mu\phi, (1 - \mu)\phi)\\ | \\\mu\\ | \\\mu(1 - \mu)/(1 + \phi)\\ |
| 2 | mean \\\mu \in (0, 1)\\, dispersion \\\phi \in (0, 1)\\ | \\(\mu\tau, (1 - \mu)\tau)\\, \\\tau = (1 - \phi)/\phi\\ | \\\mu\\ | \\\phi\\\mu(1 - \mu)\\ |

The three describe the same family: \\a + b\\ is the precision, and the
dispersion of `repar = 2` is \$\$\phi = \frac{\mathrm{Var}\[Y\]}{\mu(1 -
\mu)} = \frac{1}{1 + a + b},\$\$ the share of the largest possible
variance \\\mu(1 - \mu)\\, which is approached as \\a + b \to 0\\. It is
not a coefficient of variation.

In Lopes (2023) these are eq. `eqn_beta_p1` (shapes \\p, q\\),
"parametrizacao 1" (Ferrari and Cribari-Neto, 2004; eq. `eqn_beta_p2`)
and "parametrizacao 2" (Bayer, 2011; eq. `eqn_beta_p3`), where the
dispersion is written \\\sigma\\. The package names the second parameter
`phi` in every scheme. Regression on the shapes (`repar = 0`) is a
package extension. Admissible links per scheme: section
'Reparameterizations and links' of
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md).

## References

Lopes, J. E. (2023). *Modelos de regressao beta para dados de escala*.
Master's dissertation, Universidade Federal do Parana, Curitiba. URI:
https://hdl.handle.net/1884/86624.

Ferrari, S. L. P., and Cribari-Neto, F. (2004). Beta regression for
modelling rates and proportions. *Journal of Applied Statistics*,
**31**(7), 799–815.
[doi:10.1080/0266476042000214501](https://doi.org/10.1080/0266476042000214501)

Bayer, F. M. (2011). *Modelagem e inferencia em regressao beta*. PhD
thesis, Universidade Federal de Pernambuco.

## Examples

``` r
# One beta distribution in the three parameterisations: shapes (6, 14)
brs_repar(mu = 0.3, phi = 20, repar = 1)      # mean 0.3, precision 20
#>   shape1 shape2
#> 1      6     14
brs_repar(mu = 0.3, phi = 1 / 21, repar = 2)  # mean 0.3, dispersion 1/(1 + 20)
#>   shape1 shape2
#> 1      6     14
brs_repar(mu = 6, phi = 14, repar = 0)        # shapes p = 6, q = 14
#>   shape1 shape2
#> 1      6     14

# Back from shapes: mean, precision a + b, dispersion 1 / (1 + a + b)
sh <- brs_repar(mu = 0.3, phi = 20, repar = 1)
c(mean = sh$shape1 / (sh$shape1 + sh$shape2),
  precision = sh$shape1 + sh$shape2,
  dispersion = 1 / (1 + sh$shape1 + sh$shape2))
#>        mean   precision  dispersion 
#>  0.30000000 20.00000000  0.04761905 

# E[Y] and Var[Y] agree across parameterisations: 0.3 and 0.21 / 21 = 0.01
mu <- 0.3
c(var_repar1 = mu * (1 - mu) / (1 + 20),
  var_repar2 = mu * (1 - mu) * (1 / 21),
  var_shapes = 6 * 14 / ((6 + 14)^2 * (6 + 14 + 1)))
#> var_repar1 var_repar2 var_shapes 
#>       0.01       0.01       0.01 

# Vectorised: one row per observation
brs_repar(mu = c(0.2, 0.5, 0.8), phi = 0.1, repar = 2)
#>   shape1 shape2
#> 1    1.8    7.2
#> 2    4.5    4.5
#> 3    7.2    1.8
```
