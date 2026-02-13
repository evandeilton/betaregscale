# Reparameterize (mu, phi) into beta shape parameters

Converts a mean–dispersion pair \\(\mu, \phi)\\ to the shape parameters
\\(a, b)\\ of the beta distribution under one of three
reparameterization schemes.

## Usage

``` r
beta_reparam(mu, phi, repar = 2L)
```

## Arguments

- mu:

  Numeric vector of mean values in \\(0, 1)\\.

- phi:

  Numeric vector (or scalar) of dispersion values.

- repar:

  Integer (0, 1, or 2) selecting the scheme.

## Value

A `data.frame` with columns `shape1` and `shape2`.

## Details

The three schemes are:

- `repar = 0`:

  Direct: \\a = \mu,\\ b = \phi\\.

- `repar = 1`:

  Ferrari–Cribari-Neto: \\a = \mu\phi,\\ b = (1 - \mu)\phi\\, where
  \\\phi\\ acts as a precision parameter.

- `repar = 2`:

  Mean–variance: \\a = \mu(1-\phi)/\phi,\\ b = (1-\mu)(1-\phi)/\phi\\,
  where \\\phi \in (0,1)\\ is analogous to a coefficient of variation.

## Examples

``` r
beta_reparam(mu = 0.5, phi = 0.2, repar = 2)
#>   shape1 shape2
#> 1      2      2
```
