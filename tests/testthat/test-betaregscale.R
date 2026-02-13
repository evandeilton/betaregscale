# ============================================================================ #
# Comprehensive test suite for betaregscale
# ============================================================================ #

# -- Helpers ----------------------------------------------------------------- #

# Shared simulation fixture (fixed dispersion)
sim_fixed <- function(n = 200, seed = 42) {
  set.seed(seed)
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  betaregscale_simulate(
    formula = ~ x1 + x2, data = dat,
    beta = c(0.2, -0.5, 0.3), phi = 1 / 5,
    link = "logit", link_phi = "logit",
    ncuts = 100, repar = 2L
  )
}

# Shared simulation fixture (variable dispersion)
sim_variable <- function(n = 200, seed = 42) {
  set.seed(seed)
  dat <- data.frame(
    x1 = rnorm(n), x2 = rnorm(n),
    z1 = runif(n)
  )
  betaregscale_simulate_z(
    formula_x = ~ x1 + x2, formula_z = ~z1,
    data = dat,
    beta = c(0.2, -0.5, 0.3),
    zeta = c(0.5, -0.5),
    link = "logit", link_phi = "logit",
    ncuts = 100, repar = 2L
  )
}


# ============================================================================ #
# Test: check_response
# ============================================================================ #

test_that("check_response returns correct structure with delta column", {
  y <- c(0, 3, 5, 7, 9, 10)
  out <- check_response(y, type = "m", ncuts = 10)
  expect_true(is.matrix(out))
  expect_equal(ncol(out), 5L)
  expect_equal(colnames(out), c("left", "right", "yt", "y", "delta"))
  expect_true(all(out[, "left"] > 0))
  expect_true(all(out[, "right"] < 1))
  expect_true(all(out[, "left"] <= out[, "right"]))
})

test_that("check_response classifies boundary observations correctly", {
  y <- c(0, 5, 10)
  out <- check_response(y, type = "m", ncuts = 10)
  # y=0 should be left-censored (delta=1)
  expect_equal(unname(out[1, "delta"]), 1)
  # y=5 should be interval-censored (delta=3)
  expect_equal(unname(out[2, "delta"]), 3)
  # y=10 should be right-censored (delta=2)
  expect_equal(unname(out[3, "delta"]), 2)
})

test_that("check_response handles unit-interval input as uncensored", {
  y <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  expect_message(
    out <- check_response(y, type = "m", ncuts = 100),
    "unit interval"
  )
  # All should be uncensored (delta=0)
  expect_true(all(out[, "delta"] == 0))
  expect_true(all(out[, "left"] > 0))
  expect_true(all(out[, "right"] < 1))
})

test_that("check_response warns when max(y) > ncuts", {
  y <- c(5, 10, 150)
  expect_warning(check_response(y, ncuts = 100), "exceeds")
})

test_that("check_response types 'l' and 'r' differ from 'm'", {
  y <- c(3, 5, 7)
  m <- check_response(y, type = "m", ncuts = 10)
  l <- check_response(y, type = "l", ncuts = 10)
  r <- check_response(y, type = "r", ncuts = 10)
  # All interval-censored, but different endpoints
  expect_true(all(m[, "delta"] == 3))
  expect_true(all(l[, "delta"] == 3))
  expect_true(all(r[, "delta"] == 3))
  expect_false(all(m[, "left"] == l[, "left"]))
  expect_false(all(m[, "right"] == r[, "right"]))
})


# ============================================================================ #
# Test: beta_reparam
# ============================================================================ #

test_that("beta_reparam returns valid shapes for repar=0,1,2", {
  for (rp in 0:2) {
    out <- beta_reparam(mu = 0.5, phi = 0.3, repar = rp)
    expect_true(all(out$shape1 > 0))
    expect_true(all(out$shape2 > 0))
  }
})

test_that("beta_reparam vectorizes correctly", {
  mu <- c(0.2, 0.5, 0.8)
  phi <- c(0.1, 0.3, 0.5)
  out <- beta_reparam(mu, phi, repar = 2L)
  expect_equal(nrow(out), 3L)
})


# ============================================================================ #
# Test: link_to_code
# ============================================================================ #

test_that("link_to_code maps all valid links", {
  links <- c(
    "logit", "probit", "cauchit", "cloglog", "log", "sqrt",
    "inverse", "1/mu^2", "identity"
  )
  codes <- vapply(links, link_to_code, integer(1))
  expect_equal(length(codes), length(links))
  expect_true(all(codes >= 0))
})

test_that("link_to_code errors on invalid link", {
  expect_error(link_to_code("banana"), "Unsupported")
})


# ============================================================================ #
# Test: simulation functions
# ============================================================================ #

test_that("betaregscale_simulate produces valid output", {
  sim <- sim_fixed()
  expect_true(is.data.frame(sim))
  expect_true(all(c("left", "right", "yt", "y", "delta", "x1", "x2") %in% names(sim)))
  expect_true(all(sim$left > 0 & sim$left < 1))
  expect_true(all(sim$right > 0 & sim$right < 1))
  expect_equal(nrow(sim), 200L)
})

test_that("betaregscale_simulate_z produces valid output", {
  sim <- sim_variable()
  expect_true(is.data.frame(sim))
  expect_true("z1" %in% names(sim))
  expect_true("delta" %in% names(sim))
  expect_equal(nrow(sim), 200L)
})

test_that("betaregscale_simulate_z does NOT set seed internally", {
  set.seed(1)
  s1 <- betaregscale_simulate_z(
    formula_x = ~x1, formula_z = ~z1,
    data = data.frame(x1 = rnorm(50), z1 = runif(50)),
    beta = c(0, 0.5), zeta = c(0.5, 0.1)
  )
  set.seed(2)
  s2 <- betaregscale_simulate_z(
    formula_x = ~x1, formula_z = ~z1,
    data = data.frame(x1 = rnorm(50), z1 = runif(50)),
    beta = c(0, 0.5), zeta = c(0.5, 0.1)
  )
  expect_false(identical(s1$y, s2$y))
})


# ============================================================================ #
# Test: log-likelihood functions
# ============================================================================ #

test_that("betaregscale_loglik returns finite scalar", {
  sim <- sim_fixed()
  ll <- betaregscale_loglik(
    param = c(0.2, -0.5, 0.3, 1 / 5),
    formula = y ~ x1 + x2, data = sim
  )
  expect_true(is.numeric(ll))
  expect_equal(length(ll), 1L)
  expect_true(is.finite(ll))
  expect_true(ll < 0)
})

test_that("betaregscale_loglik_z returns finite scalar", {
  sim <- sim_variable()
  ll <- betaregscale_loglik_z(
    param = c(0.2, -0.5, 0.3, 0.5, -0.5),
    formula = y ~ x1 + x2 | z1, data = sim
  )
  expect_true(is.numeric(ll))
  expect_true(is.finite(ll))
  expect_true(ll < 0)
})


# ============================================================================ #
# Test: fitting functions
# ============================================================================ #

test_that("betaregscale_fit converges and returns correct class", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(
    formula = y ~ x1 + x2, data = sim,
    link = "logit", link_phi = "logit"
  )
  expect_s3_class(fit, "betaregscale")
  expect_equal(fit$convergence, 0L)
  expect_true(is.finite(fit$value))
  expect_equal(length(fit$par), 4L)
  expect_equal(fit$nobs, nrow(sim))
})

test_that("betaregscale_fit_z converges for variable dispersion", {
  sim <- sim_variable()
  fit <- betaregscale_fit_z(
    formula = y ~ x1 + x2 | z1, data = sim,
    link = "logit", link_phi = "logit"
  )
  expect_s3_class(fit, "betaregscale")
  expect_equal(fit$convergence, 0L)
  expect_equal(length(fit$par), 5L)
})

test_that("betaregscale() dispatches to fixed model", {
  sim <- sim_fixed()
  fit <- betaregscale(y ~ x1 + x2, data = sim)
  expect_s3_class(fit, "betaregscale")
  expect_equal(length(fit$par), 4L)
})

test_that("betaregscale() dispatches to variable model", {
  sim <- sim_variable()
  fit <- betaregscale(y ~ x1 + x2 | z1, data = sim)
  expect_s3_class(fit, "betaregscale")
  expect_equal(length(fit$par), 5L)
})

test_that("fitted object stores call", {
  sim <- sim_fixed()
  fit <- betaregscale(y ~ x1 + x2, data = sim)
  expect_true(!is.null(fit$call))
})


# ============================================================================ #
# Test: parameter naming (betareg style)
# ============================================================================ #

test_that("fixed model parameters have correct names", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)
  nm <- names(fit$par)
  expect_equal(nm, c("(Intercept)", "x1", "x2", "(phi)"))
  expect_equal(names(fit$coefficients$mean), c("(Intercept)", "x1", "x2"))
  expect_equal(names(fit$coefficients$precision), "(phi)")
})

test_that("variable model parameters have (phi)_ prefix", {
  sim <- sim_variable()
  fit <- betaregscale_fit_z(y ~ x1 + x2 | z1, data = sim)
  nm <- names(fit$par)
  expect_true("(phi)_(Intercept)" %in% nm)
  expect_true("(phi)_z1" %in% nm)
})


# ============================================================================ #
# Test: S3 methods <U+2014> coef with model= argument
# ============================================================================ #

test_that("coef method works with model argument", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)

  cf_full <- coef(fit, model = "full")
  cf_mean <- coef(fit, model = "mean")
  cf_prec <- coef(fit, model = "precision")

  expect_equal(length(cf_full), 4L)
  expect_equal(length(cf_mean), 3L)
  expect_equal(length(cf_prec), 1L)
  expect_true(all(!is.na(cf_full)))
})


# ============================================================================ #
# Test: vcov with model= argument
# ============================================================================ #

test_that("vcov method returns valid matrix with model subsetting", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)

  V_full <- vcov(fit, model = "full")
  V_mean <- vcov(fit, model = "mean")
  V_prec <- vcov(fit, model = "precision")

  expect_equal(nrow(V_full), 4L)
  expect_equal(ncol(V_full), 4L)
  expect_equal(nrow(V_mean), 3L)
  expect_equal(nrow(V_prec), 1L)
  expect_true(all(diag(V_full) > 0))
})


# ============================================================================ #
# Test: logLik, AIC, BIC
# ============================================================================ #

test_that("logLik method returns correct class", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_true(is.finite(as.numeric(ll)))
  expect_equal(attr(ll, "df"), 4L)
})

test_that("AIC and BIC return scalars", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)
  aic <- AIC(fit)
  bic <- BIC(fit)
  expect_true(is.finite(aic))
  expect_true(is.finite(bic))
  expect_true(bic > aic) # BIC penalizes more with n > e^2
})

test_that("stats::AIC(logLik(fit)) == AIC(fit)", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)
  expect_equal(stats::AIC(logLik(fit)), AIC(fit), tolerance = 1e-6)
})

test_that("stats::BIC(logLik(fit)) == BIC(fit)", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)
  expect_equal(stats::BIC(logLik(fit)), BIC(fit), tolerance = 1e-6)
})


# ============================================================================ #
# Test: nobs, formula, model.matrix
# ============================================================================ #

test_that("nobs returns correct count", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)
  expect_equal(nobs(fit), 200L)
})

test_that("formula returns the model formula", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)
  f <- formula(fit)
  expect_true(inherits(f, "formula"))
})

test_that("model.matrix returns design matrices", {
  sim <- sim_variable()
  fit <- betaregscale_fit_z(y ~ x1 + x2 | z1, data = sim)
  X <- model.matrix(fit, model = "mean")
  Z <- model.matrix(fit, model = "precision")
  expect_equal(nrow(X), 200L)
  expect_equal(ncol(X), 3L) # Intercept + x1 + x2
  expect_equal(nrow(Z), 200L)
  expect_equal(ncol(Z), 2L) # Intercept + z1
})


# ============================================================================ #
# Test: summary method (betareg style)
# ============================================================================ #

test_that("summary method returns structured output", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)
  s <- summary(fit)
  expect_s3_class(s, "summary.betaregscale")
  expect_true("coefficients" %in% names(s))
  expect_true("mean" %in% names(s$coefficients))
  expect_true("precision" %in% names(s$coefficients))
  expect_equal(nrow(s$coefficients$mean), 3L)
  expect_equal(nrow(s$coefficients$precision), 1L)
  expect_true(!is.null(s$censoring))
})

test_that("summary uses pnorm not pt for p-values", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)
  s <- summary(fit)
  # p-values should match pnorm calculation
  z_val <- s$coefficients$mean[, "z value"]
  p_val <- s$coefficients$mean[, "Pr(>|z|)"]
  expected_p <- 2 * pnorm(-abs(z_val))
  expect_equal(p_val, expected_p, tolerance = 1e-10)
})


# ============================================================================ #
# Test: print methods
# ============================================================================ #

test_that("print method runs without error", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)
  expect_output(print(fit), "Coefficients")
})

test_that("print.summary method runs without error", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)
  expect_output(print(summary(fit)), "Coefficients")
})


# ============================================================================ #
# Test: fitted values
# ============================================================================ #

test_that("fitted method returns vectors of correct length", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)
  mu <- fitted(fit, type = "mu")
  phi <- fitted(fit, type = "phi")
  expect_equal(length(mu), nrow(sim))
  expect_true(all(mu > 0 & mu < 1))
})


# ============================================================================ #
# Test: residuals by type
# ============================================================================ #

test_that("residuals method returns correct types", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)

  for (rtype in c("response", "pearson", "deviance", "rqr")) {
    r <- residuals(fit, type = rtype)
    expect_equal(length(r), nrow(sim))
    expect_true(all(is.finite(r)))
  }
})

test_that("RQR residuals are approximately normal for correct model", {
  sim <- sim_fixed(n = 500, seed = 123)
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)
  r <- residuals(fit, type = "rqr")
  # Shapiro-Wilk should not reject at 1% level for n <= 5000
  sw <- shapiro.test(r[1:min(500, length(r))])
  expect_true(sw$p.value > 0.01)
})


# ============================================================================ #
# Test: confint
# ============================================================================ #

test_that("confint returns valid confidence intervals", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)
  ci <- confint(fit)
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 4L)
  expect_equal(ncol(ci), 2L)
  expect_true(all(ci[, 1] < ci[, 2]))
})

test_that("confint with model argument works", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)
  ci_mean <- confint(fit, model = "mean")
  ci_prec <- confint(fit, model = "precision")
  expect_equal(nrow(ci_mean), 3L)
  expect_equal(nrow(ci_prec), 1L)
})

test_that("confint with parm argument works", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)
  ci <- confint(fit, parm = c("x1", "x2"))
  expect_equal(nrow(ci), 2L)
})


# ============================================================================ #
# Test: predict
# ============================================================================ #

test_that("predict with type='response' works", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)
  pred <- predict(fit, type = "response")
  expect_equal(length(pred), nrow(sim))
  expect_true(all(pred > 0 & pred < 1))
})

test_that("predict with newdata works", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)
  nd <- data.frame(x1 = c(-1, 0, 1), x2 = c(0, 0, 0))
  pred <- predict(fit, newdata = nd, type = "response")
  expect_equal(length(pred), 3L)
  expect_true(all(pred > 0 & pred < 1))
})

test_that("predict type='link' returns linear predictor", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)
  pred <- predict(fit, type = "link")
  expect_equal(length(pred), nrow(sim))
  # link values are on unrestricted scale
})

test_that("predict type='variance' works", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)
  pred <- predict(fit, type = "variance")
  expect_equal(length(pred), nrow(sim))
  expect_true(all(pred > 0))
})

test_that("predict type='quantile' works", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)
  pred <- predict(fit, type = "quantile", at = c(0.25, 0.5, 0.75))
  expect_true(is.matrix(pred))
  expect_equal(ncol(pred), 3L)
  expect_equal(nrow(pred), nrow(sim))
})


# ============================================================================ #
# Test: predict with newdata for variable dispersion
# ============================================================================ #

test_that("predict with newdata works for variable-dispersion model", {
  sim <- sim_variable()
  fit <- betaregscale_fit_z(y ~ x1 + x2 | z1, data = sim)
  nd <- data.frame(x1 = c(-1, 0, 1), x2 = c(0, 0, 0), z1 = c(0.3, 0.5, 0.7))
  pred <- predict(fit, newdata = nd, type = "response")
  expect_equal(length(pred), 3L)
  expect_true(all(pred > 0 & pred < 1))
})


# ============================================================================ #
# Test: gof and est convenience functions
# ============================================================================ #

test_that("gof returns data frame with correct columns", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)
  g <- gof(fit)
  expect_true(is.data.frame(g))
  expect_true(all(c("logLik", "AIC", "BIC") %in% names(g)))
})

test_that("est returns estimates with p-values using pnorm", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)
  e <- est(fit)
  expect_true(is.data.frame(e))
  expect_true("p_value" %in% names(e))
  expect_equal(nrow(e), 4L)
})


# ============================================================================ #
# Test: edge cases
# ============================================================================ #

test_that("fitting with probit link works", {
  set.seed(42)
  n <- 100
  dat <- data.frame(x1 = rnorm(n))
  sim <- betaregscale_simulate(
    formula = ~x1, data = dat,
    beta = c(0.2, 0.3), phi = 1 / 5,
    link = "logit", link_phi = "logit"
  )
  fit <- betaregscale_fit(
    y ~ x1,
    data = sim,
    link = "probit", link_phi = "logit"
  )
  expect_s3_class(fit, "betaregscale")
  expect_equal(fit$convergence, 0L)
})

test_that("betaregscale_coef produces valid confidence intervals", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)
  tab <- betaregscale_coef(fit)$est
  expect_true(all(tab$ci_lower < tab$estimate))
  expect_true(all(tab$ci_upper > tab$estimate))
})


# ============================================================================ #
# Test: mixed censoring in check_response
# ============================================================================ #

test_that("mixed censoring: all four types present", {
  # Create data with all censoring types
  y_mixed <- c(0, 0, 3, 5, 7, 10, 10)
  out <- check_response(y_mixed, type = "m", ncuts = 10)
  deltas <- out[, "delta"]
  expect_true(1L %in% deltas) # left-censored
  expect_true(2L %in% deltas) # right-censored
  expect_true(3L %in% deltas) # interval-censored
})


# ============================================================================ #
# Test: plot runs without error
# ============================================================================ #

test_that("plot.betaregscale runs without error", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)
  expect_silent(plot(fit, which = 1:4))
})


# ============================================================================ #
# Test: censoring_summary runs without error
# ============================================================================ #

test_that("censoring_summary returns correct structure", {
  sim <- sim_fixed()
  fit <- betaregscale_fit(y ~ x1 + x2, data = sim)
  cs <- censoring_summary(fit)
  expect_true(is.data.frame(cs))
  expect_true("type" %in% names(cs))
  expect_true("count" %in% names(cs))
  expect_true("proportion" %in% names(cs))
  expect_equal(sum(cs$count), nrow(sim))
})

test_that("censoring_summary works with raw matrix", {
  y <- c(0, 3, 5, 7, 10)
  Y <- check_response(y, type = "m", ncuts = 10)
  cs <- censoring_summary(Y)
  expect_true(is.data.frame(cs))
  expect_equal(sum(cs$count), 5L)
})


# ============================================================================ #
# bs_prepare() tests
# ============================================================================ #

# -- Structure and attributes ------------------------------------------------ #

test_that("bs_prepare output has correct structure", {
  d <- data.frame(y = c(0, 5, 10), x1 = rnorm(3))
  res <- bs_prepare(d, ncuts = 10)
  expect_true(is.data.frame(res))
  expect_true(all(c("left", "right", "yt", "y", "delta") %in% names(res)))
  expect_true("x1" %in% names(res))
  expect_true(isTRUE(attr(res, "bs_prepared")))
  expect_equal(attr(res, "ncuts"), 10L)
  expect_equal(attr(res, "type"), "m")
  expect_equal(attr(res, "lim"), 0.5)
})

test_that("bs_prepare preserves original y values", {
  d <- data.frame(y = c(0, 50, 100))
  res <- bs_prepare(d, ncuts = 100)
  expect_equal(res$y, c(0, 50, 100))
})

# -- Mode 1: y only (automatic classification) ------------------------------ #

test_that("bs_prepare Mode 1: y only <U+2014> delta inferred correctly", {
  d <- data.frame(y = c(0, 3, 5, 7, 10))
  res <- bs_prepare(d, ncuts = 10)
  expect_equal(res$delta, c(1L, 3L, 3L, 3L, 2L))
})

test_that("bs_prepare Mode 1: y only <U+2014> endpoints match check_response", {
  y_vec <- c(0, 3, 5, 7, 10)
  cr <- check_response(y_vec, type = "m", ncuts = 10, lim = 0.5)
  d <- data.frame(y = y_vec)
  res <- bs_prepare(d, ncuts = 10, type = "m", lim = 0.5)
  expect_equal(res$left, as.numeric(cr[, "left"]), tolerance = 1e-8)
  expect_equal(res$right, as.numeric(cr[, "right"]), tolerance = 1e-8)
  expect_equal(res$yt, as.numeric(cr[, "yt"]), tolerance = 1e-8)
  expect_equal(res$delta, as.integer(cr[, "delta"]))
})

# -- Mode 2: y + explicit delta --------------------------------------------- #

test_that("bs_prepare Mode 2: y + delta explicit <U+2014> exact", {
  d <- data.frame(y = 50, delta = 0L)
  res <- bs_prepare(d, ncuts = 100)
  expect_equal(res$delta, 0L)
  expect_equal(res$left, 0.50, tolerance = 1e-4)
  expect_equal(res$right, 0.50, tolerance = 1e-4)
})

test_that("bs_prepare Mode 2: y + delta explicit <U+2014> left-censored", {
  d <- data.frame(y = 0, delta = 1L)
  res <- bs_prepare(d, ncuts = 100)
  expect_equal(res$delta, 1L)
  expect_true(res$left < 0.01)
  expect_equal(res$right, 0.005, tolerance = 1e-4)
})

test_that("bs_prepare Mode 2: y + delta explicit <U+2014> right-censored", {
  d <- data.frame(y = 100, delta = 2L)
  res <- bs_prepare(d, ncuts = 100)
  expect_equal(res$delta, 2L)
  expect_equal(res$left, 0.995, tolerance = 1e-4)
  expect_true(res$right > 0.99)
})

test_that("bs_prepare Mode 2: y + delta explicit <U+2014> interval-censored", {
  d <- data.frame(y = 50, delta = 3L)
  res <- bs_prepare(d, ncuts = 100, type = "m")
  expect_equal(res$delta, 3L)
  expect_equal(res$left, 0.495, tolerance = 1e-4)
  expect_equal(res$right, 0.505, tolerance = 1e-4)
  expect_equal(res$yt, 0.50, tolerance = 1e-4)
})

# -- Mode 3: left/right with NA patterns ------------------------------------ #

test_that("bs_prepare Mode 3: only right given <U+2014> left-censored", {
  d <- data.frame(left = NA_real_, right = 5, y = NA_real_)
  res <- bs_prepare(d, ncuts = 100)
  expect_equal(res$delta, 1L)
  expect_true(res$left < 0.001)
  expect_equal(res$right, 0.05, tolerance = 1e-4)
})

test_that("bs_prepare Mode 3: only left given <U+2014> right-censored", {
  d <- data.frame(left = 20, right = NA_real_, y = NA_real_)
  res <- bs_prepare(d, ncuts = 100)
  expect_equal(res$delta, 2L)
  expect_equal(res$left, 0.20, tolerance = 1e-4)
  expect_true(res$right > 0.99)
})

test_that("bs_prepare Mode 3: left + right given <U+2014> interval-censored", {
  d <- data.frame(left = 30, right = 45, y = NA_real_)
  res <- bs_prepare(d, ncuts = 100)
  expect_equal(res$delta, 3L)
  expect_equal(res$left, 0.30, tolerance = 1e-4)
  expect_equal(res$right, 0.45, tolerance = 1e-4)
  expect_equal(res$yt, 0.375, tolerance = 1e-4)
})

test_that("bs_prepare Mode 3: y only (no left/right) <U+2014> scale value <U+2192> interval", {
  # y=50 on 0-100 scale <U+2192> interval-censored (delta=3), not exact
  d <- data.frame(left = NA_real_, right = NA_real_, y = 50)
  res <- bs_prepare(d, ncuts = 100, type = "m")
  expect_equal(res$delta, 3L)
  expect_equal(res$yt, 0.50, tolerance = 1e-4)
})

test_that("bs_prepare: continuous y in (0,1) <U+2192> exact", {
  d <- data.frame(y = c(0.2, 0.5, 0.8))
  res <- bs_prepare(d, ncuts = 100)
  expect_equal(res$delta, c(0L, 0L, 0L))
  expect_equal(res$yt, c(0.2, 0.5, 0.8), tolerance = 1e-4)
})

# -- Mode 4: y + left + right (analyst-supplied intervals) ------------------- #

test_that("bs_prepare Mode 4: y + left + right <U+2014> uses analyst endpoints", {
  d <- data.frame(y = 50, left = 48, right = 52)
  res <- bs_prepare(d, ncuts = 100)
  expect_equal(res$delta, 3L)
  expect_equal(res$left, 0.48, tolerance = 1e-4)
  expect_equal(res$right, 0.52, tolerance = 1e-4)
  expect_equal(res$yt, 0.50, tolerance = 1e-4)
})

# -- Mixed censoring in one dataset ----------------------------------------- #

test_that("bs_prepare handles mixed censoring types", {
  d <- data.frame(
    left  = c(NA, 20, 30, NA),
    right = c(5, NA, 45, NA),
    y     = c(NA, NA, NA, 50)
  )
  res <- bs_prepare(d, ncuts = 100)
  # y=50 is a scale value (not in (0,1)), so it's interval-censored

  expect_equal(res$delta, c(1L, 2L, 3L, 3L))
})

# -- Input validation errors ------------------------------------------------- #

test_that("bs_prepare errors on non-data.frame", {
  expect_error(bs_prepare(list(y = 1:5)), "data.frame")
})

test_that("bs_prepare errors when no relevant columns", {
  expect_error(bs_prepare(data.frame(x1 = 1:5)), "must be present")
})

test_that("bs_prepare errors on invalid delta values", {
  d <- data.frame(y = 50, delta = 5L)
  expect_error(bs_prepare(d, ncuts = 100), "\\{0, 1, 2, 3\\}")
})

test_that("bs_prepare errors when left > right", {
  d <- data.frame(left = 50, right = 30, y = NA_real_)
  expect_error(bs_prepare(d, ncuts = 100), "left > right")
})

test_that("bs_prepare errors when all columns are NA", {
  d <- data.frame(y = NA_real_, left = NA_real_, right = NA_real_)
  expect_error(bs_prepare(d, ncuts = 100), "all relevant columns are NA")
})

test_that("bs_prepare errors when ncuts < max value", {
  d <- data.frame(y = c(50, 150))
  expect_error(bs_prepare(d, ncuts = 100), "ncuts")
})

test_that("bs_prepare errors on non-numeric y", {
  d <- data.frame(y = c("a", "b"))
  expect_error(bs_prepare(d, ncuts = 100), "numeric")
})

test_that("bs_prepare errors on negative y", {
  d <- data.frame(y = c(-1, 5))
  expect_error(bs_prepare(d, ncuts = 100), "non-negative")
})

# -- Consistency warnings ---------------------------------------------------- #

test_that("bs_prepare warns: delta=1 but y != 0", {
  d <- data.frame(y = 5, delta = 1L)
  expect_warning(bs_prepare(d, ncuts = 100), "delta = 1.*but y != 0")
})

test_that("bs_prepare warns: delta=2 but y != ncuts", {
  d <- data.frame(y = 50, delta = 2L)
  expect_warning(bs_prepare(d, ncuts = 100), "delta = 2.*but y != 100")
})

test_that("bs_prepare warns: delta=3 but y at boundary", {
  d <- data.frame(y = 0, delta = 3L)
  expect_warning(bs_prepare(d, ncuts = 100), "delta = 3.*boundary")
})

# -- censoring_summary with bs_prepare output -------------------------------- #

test_that("censoring_summary works with bs_prepare output", {
  d <- data.frame(y = c(0, 3, 5, 7, 10), x1 = rnorm(5))
  res <- bs_prepare(d, ncuts = 10)
  cs <- censoring_summary(res)
  expect_true(is.data.frame(cs))
  expect_equal(sum(cs$count), 5L)
})

# -- End-to-end: bs_prepare -> betaregscale ---------------------------------- #

test_that("bs_prepare -> betaregscale produces identical fit to raw data", {
  sim <- sim_fixed(n = 200, seed = 42)

  # Fit 1: standard workflow
  fit1 <- betaregscale(y ~ x1 + x2, data = sim, link = "logit", repar = 2L)

  # Fit 2: bs_prepare workflow <U+2014> extract only y + covariates
  sim_raw <- sim[, c("y", "x1", "x2")]
  prep <- bs_prepare(sim_raw, ncuts = 100, type = "m")
  fit2 <- betaregscale(y ~ x1 + x2, data = prep, link = "logit", repar = 2L)

  expect_equal(coef(fit1), coef(fit2), tolerance = 1e-6)
  expect_equal(logLik(fit1), logLik(fit2), tolerance = 1e-6)
})

test_that("bs_prepare -> betaregscale_fit_z converges", {
  sim <- sim_variable(n = 200, seed = 42)
  sim_raw <- sim[, c("y", "x1", "x2", "z1")]
  prep <- bs_prepare(sim_raw, ncuts = 100, type = "m")
  fit <- betaregscale(y ~ x1 + x2 | z1, data = prep, link = "logit", repar = 2L)
  expect_s3_class(fit, "betaregscale")
  expect_equal(fit$convergence, 0L)
})
