# Wald confidence intervals

Wald intervals \\\hat\theta_j \pm z\_{1 - \alpha/2} SE_j\\ on the link
scale, with \\SE_j\\ from
[`vcov.brs`](https://evandeilton.github.io/betaregscale/reference/vcov.brs.md)
(Lopes, 2023, "Inferencia").

## Usage

``` r
# S3 method for class 'brs'
confint(
  object,
  parm,
  level = 0.95,
  model = c("full", "mean", "precision"),
  ...
)
```

## Arguments

- object:

  A fitted `"brs"` object.

- parm:

  Character or integer: which parameters. If missing, all parameters of
  `model` are returned.

- level:

  Confidence level (default 0.95).

- model:

  Character: `"full"`, `"mean"` or `"precision"`.

- ...:

  Currently ignored.

## Value

Matrix with the lower and upper limits.

## Details

Intervals for a mean or precision on the response scale follow by the
inverse link of the limits (monotone links). A limit is `NA` when the
variance is not estimable (see 'Fit diagnostics' in
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md)).
For small samples or parameters near the border of the scale,
[`brs_bootstrap`](https://evandeilton.github.io/betaregscale/reference/brs_bootstrap.md)
gives intervals that do not rely on the normal approximation.

## References

Lopes, J. E. (2023). *Modelos de regressao beta para dados de escala*.
Master's dissertation, Universidade Federal do Parana, Curitiba. URI:
https://hdl.handle.net/1884/86624.

## See also

[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md),
[`vcov.brs`](https://evandeilton.github.io/betaregscale/reference/vcov.brs.md),
[`brs_bootstrap`](https://evandeilton.github.io/betaregscale/reference/brs_bootstrap.md),
[`brs_est`](https://evandeilton.github.io/betaregscale/reference/brs_est.md)

## Examples

``` r
set.seed(2023)
d <- data.frame(x = runif(150))
s <- brs_sim(~ x, data = d, beta = c(-0.5, 1), phi = qlogis(0.3), ncuts = 10)
fit <- brs(y ~ x, data = s)
confint(fit)
#>                  2.5 %     97.5 %
#> (Intercept) -0.5662544  0.1735143
#> x           -0.2500579  0.9946682
#> (phi)       -1.0194170 -0.5981087
# Mean at x = 0 on (0, 1): inverse logit of the intercept limits
plogis(confint(fit, parm = "(Intercept)"))
#>                 2.5 %    97.5 %
#> (Intercept) 0.3621016 0.5432701
```
