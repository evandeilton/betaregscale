# Fit a variable-dispersion beta interval regression model

Estimates the parameters of a beta regression model with
observation-specific dispersion governed by a second linear predictor.
Both submodels are estimated jointly via maximum likelihood, using the
complete likelihood with mixed censoring (Lopes, 2024, Eq. 2.24).

## Usage

``` r
brs_fit_var(
  formula,
  data,
  link = "logit",
  link_phi = "logit",
  hessian_method = c("numDeriv", "optim"),
  ncuts = 100L,
  lim = 0.5,
  repar = 2L,
  method = c("BFGS", "L-BFGS-B")
)
```

## Arguments

- formula:

  A [`Formula`](https://rdrr.io/pkg/Formula/man/Formula.html)-style
  formula with two parts: `y ~ x1 + x2 | z1 + z2`.

- data:

  Data frame.

- link:

  Mean link function (default `"logit"`).

- link_phi:

  Dispersion link function (default `"logit"`).

- hessian_method:

  Character: `"numDeriv"` or `"optim"`.

- ncuts:

  Number of scale categories (default 100).

- lim:

  Uncertainty half-width (default 0.5).

- repar:

  Reparameterization scheme (default 2).

- method:

  Optimization method (default `"BFGS"`).

## Value

An object of class `"brs"`.

## References

Hawker, G. A., Mian, S., Kendzerska, T., and French, M. (2011). Measures
of adult pain: Visual Analog Scale for Pain (VAS Pain), Numeric Rating
Scale for Pain (NRS Pain), McGill Pain Questionnaire (MPQ), Short-Form
McGill Pain Questionnaire (SF-MPQ), Chronic Pain Grade Scale (CPGS),
Short Form-36 Bodily Pain Scale (SF-36 BPS), and Measure of Intermittent
and Constant Osteoarthritis Pain (ICOAP). Arthritis Care and Research,
63(S11), S240-S252. doi:10.1002/acr.20543.

Hjermstad, M. J., Fayers, P. M., Haugen, D. F., et al. (2011). Studies
comparing Numerical Rating Scales, Verbal Rating Scales, and Visual
Analogue Scales for assessment of pain intensity in adults: a systematic
literature review. Journal of Pain and Symptom Management, 41(6),
1073-1093. doi:10.1016/j.jpainsymman.2010.08.016.

## Examples

``` r
set.seed(42)
n <- 100
dat <- data.frame(
  x1 = rnorm(n), x2 = rnorm(n),
  z1 = runif(n)
)
sim <- brs_sim(
  formula = ~ x1 + x2 | z1,
  data = dat,
  beta = c(0.2, -0.5, 0.3),
  zeta = c(1, 1.2)
)
fit <- brs_fit_var(
  formula = y ~ x1 + x2 | z1, data = sim,
  link = "logit", link_phi = "logit"
)
print(fit)
#> 
#> Call:
#> brs_fit_var(formula = y ~ x1 + x2 | z1, data = sim, link = "logit", 
#>     link_phi = "logit")
#> 
#> Coefficients (mean model with logit link):
#> (Intercept)          x1          x2 
#>      0.0959     -0.3912      0.2965 
#> 
#> Phi coefficients (precision model with logit link):
#> (phi)_(Intercept)          (phi)_z1 
#>            1.4876            0.0040 
#> 
```
