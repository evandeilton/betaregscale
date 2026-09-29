# Link functions allowed under each reparameterization

Each parameter of the beta distribution is modelled through a link whose
inverse maps the real line onto the parameter's domain. Parameters on
\\(0, 1)\\ (the mean under `repar = 1, 2`; the dispersion under
`repar = 2`) use `"logit"`, `"probit"`, `"cauchit"` or `"cloglog"`.
Parameters on \\(0, \infty)\\ (both shapes under `repar = 0`; the
precision under `repar = 1`) use `"log"` or `"sqrt"`. `"identity"`,
`"inverse"` and `"1/mu^2"` are not accepted for positive parameters:
their inverse does not map the real line onto \\(0, \infty)\\ (negative
values, a discontinuity at 0, undefined for \\\eta \le 0\\).

## Usage

``` r
.brs_link_table
```
