# Fit a mixed-effects beta interval regression model

Fits a beta interval-censored mixed model with Gaussian random
intercepts/slopes using marginal maximum likelihood. The implementation
supports random-effects formulas such as `~ 1 | group` and
`~ 1 + x | group`, and offers three integration methods for the random
effects: Laplace approximation, Adaptive Gauss-Hermite Quadrature
(AGHQ), and Quasi-Monte Carlo (QMC).

## Usage

``` r
brsmm(
  formula,
  random = ~1 | id,
  data,
  link = "logit",
  link_phi = "logit",
  repar = 2L,
  ncuts = 100L,
  lim = 0.5,
  int_method = c("laplace", "aghq", "qmc"),
  n_points = 11L,
  qmc_points = 1024L,
  start = NULL,
  method = c("BFGS", "L-BFGS-B"),
  hessian_method = c("numDeriv", "optim"),
  control = list(maxit = 2000L)
)
```

## Arguments

- formula:

  Model formula. Supports one- or two-part formulas: `y ~ x1 + x2` or
  `y ~ x1 + x2 | z1 + z2`.

- random:

  Random-effects specification of the form `~ terms | group`, e.g.
  `~ 1 | id` or `~ 1 + x | id`.

- data:

  Data frame.

- link:

  Mean link function.

- link_phi:

  Precision link function.

- repar:

  Beta reparameterization code (0, 1, 2).

- ncuts:

  Number of categories on the original scale.

- lim:

  Half-width used to construct interval endpoints.

- int_method:

  Integration method: `"laplace"` (default), `"aghq"`, or `"qmc"`.

- n_points:

  Number of quadrature points for `int_method="aghq"`. Ignored for other
  methods. Default is 11.

- qmc_points:

  Number of QMC points for `int_method="qmc"`. Default is 1024.

- start:

  Optional numeric vector of starting values (`beta`, `gamma`, and
  packed lower-Cholesky random parameters).

- method:

  Optimizer passed to [`optim`](https://rdrr.io/r/stats/optim.html).

- hessian_method:

  `"numDeriv"` (default) or `"optim"`.

- control:

  Control list for [`optim`](https://rdrr.io/r/stats/optim.html).

## Value

An object of class `"brsmm"`.

## Details

The conditional contribution for each observation follows the same mixed
censoring likelihood used by
[`brs`](https://evandeilton.github.io/betaregscale/reference/brs.md):

1.  \\\delta=0\\: exact contribution via beta density,

2.  \\\delta=1\\: left-censored contribution via beta CDF,

3.  \\\delta=2\\: right-censored contribution via survival CDF,

4.  \\\delta=3\\: interval contribution via CDF difference.

For group \\i\\, the random-effects vector \\\mathbf{b}\_i \sim
N(\mathbf{0}, D)\\ is integrated out numerically.

- `"laplace"`: Uses a second-order Laplace approximation at the
  conditional mode. Fast and generally accurate for \\n_i\\ large.

- `"aghq"`: Adaptive Gauss-Hermite Quadrature. Uses `n_points`
  quadrature nodes centered and scaled by the conditional mode and
  curvature. More accurate than Laplace, especially for small \\n_i\\.

- `"qmc"`: Quasi-Monte Carlo integration using a Halton sequence. Uses
  `qmc_points` evaluation points. Suitable for high-dimensional
  integration (future proofing) or checking robustness.

## References

Ferrari, S. and Cribari-Neto, F. (2004). Beta regression for modelling
rates and proportions. *Journal of Applied Statistics*, **31**(7),
799–815.

## Examples

``` r
# \donttest{
set.seed(123)
g <- 15
ni <- 8
id <- factor(rep(seq_len(g), each = ni))
n <- length(id)
x1 <- rnorm(n)
b <- rnorm(g, sd = 0.5)
eta_mu <- 0.2 + 0.6 * x1 + b[as.integer(id)]
mu <- plogis(eta_mu)
phi <- plogis(-0.2 + 0.2 * x1)
shp <- brs_repar(mu = mu, phi = phi, repar = 2)
y <- round(stats::rbeta(n, shp$shape1, shp$shape2) * 100)
d <- data.frame(y = y, x1 = x1, id = id)

fit_mm <- brsmm(y ~ x1, random = ~ 1 | id, data = d, repar = 2)
fit_mm
#> 
#> Call:
#> brsmm(formula = y ~ x1, random = ~1 | id, data = d, repar = 2)
#> 
#> Mixed beta interval model (Laplace)
#> Observations: 120  | Groups: 15 
#> Log-likelihood: -499.3888 
#> Random SD: 0.2207 
#> Convergence code: 0 
# }
```
