# ============================================================================ #
# Comprehensive test suite for betaregscale
# ============================================================================ #

# -- Helpers ----------------------------------------------------------------- #

# Shared simulation fixture (fixed dispersion)
sim_fixed <- function(n = 200, seed = 42) {
  set.seed(seed)
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  brs_sim(
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
  brs_sim(
    formula = ~ x1 + x2 | z1,
    data = dat,
    beta = c(0.2, -0.5, 0.3),
    zeta = c(0.5, -0.5),
    link = "logit", link_phi = "logit",
    ncuts = 100, repar = 2L
  )
}


# ============================================================================ #
# Test: brs_check
# ============================================================================ #

test_that("brs_check returns correct structure with delta column", {
  y <- c(0, 3, 5, 7, 9, 10)
  out <- brs_check(y, ncuts = 10)
  expect_true(is.matrix(out))
  expect_equal(ncol(out), 5L)
  expect_equal(colnames(out), c("left", "right", "yt", "y", "delta"))
  expect_true(all(out[, "left"] > 0))
  expect_true(all(out[, "right"] < 1))
  expect_true(all(out[, "left"] <= out[, "right"]))
})

test_that("brs_check classifies boundary observations correctly", {
  y <- c(0, 5, 10)
  out <- brs_check(y, ncuts = 10)
  # y=0 should be left-censored (delta=1)
  expect_equal(unname(out[1, "delta"]), 1)
  # y=5 should be interval-censored (delta=3)
  expect_equal(unname(out[2, "delta"]), 3)
  # y=10 should be right-censored (delta=2)
  expect_equal(unname(out[3, "delta"]), 2)
})

test_that("brs_check handles unit-interval input as uncensored", {
  y <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  expect_message(
    out <- brs_check(y, ncuts = 100),
    "unit interval"
  )
  # All should be uncensored (delta=0)
  expect_true(all(out[, "delta"] == 0))
  expect_true(all(out[, "left"] > 0))
  expect_true(all(out[, "right"] < 1))
})

test_that("brs_check warns when max(y) > ncuts", {
  y <- c(5, 10, 150)
  expect_warning(brs_check(y, ncuts = 100), "exceeds")
})


# ============================================================================ #
# Test: brs_repar
# ============================================================================ #

test_that("brs_repar returns valid shapes for repar=0,1,2", {
  for (rp in 0:2) {
    out <- brs_repar(mu = 0.5, phi = 0.3, repar = rp)
    expect_true(all(out$shape1 > 0))
    expect_true(all(out$shape2 > 0))
  }
})

test_that("brs_repar vectorizes correctly", {
  mu <- c(0.2, 0.5, 0.8)
  phi <- c(0.1, 0.3, 0.5)
  out <- brs_repar(mu, phi, repar = 2L)
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

test_that("brs_sim produces valid output", {
  sim <- sim_fixed()
  expect_true(is.data.frame(sim))
  expect_true(all(c("left", "right", "yt", "y", "delta", "x1", "x2") %in% names(sim)))
  expect_true(all(sim$left > 0 & sim$left < 1))
  expect_true(all(sim$right > 0 & sim$right < 1))
  expect_equal(nrow(sim), 200L)
})

test_that("brs_sim supports variable dispersion via two-part formula", {
  sim <- sim_variable()
  expect_true(is.data.frame(sim))
  expect_true("z1" %in% names(sim))
  expect_true("delta" %in% names(sim))
  expect_equal(nrow(sim), 200L)
})

test_that("brs_sim does NOT set seed internally for variable dispersion", {
  set.seed(1)
  s1 <- brs_sim(
    formula = ~ x1 | z1,
    data = data.frame(x1 = rnorm(50), z1 = runif(50)),
    beta = c(0, 0.5), zeta = c(0.5, 0.1)
  )
  set.seed(2)
  s2 <- brs_sim(
    formula = ~ x1 | z1,
    data = data.frame(x1 = rnorm(50), z1 = runif(50)),
    beta = c(0, 0.5), zeta = c(0.5, 0.1)
  )
  expect_false(identical(s1$y, s2$y))
})


# ============================================================================ #
# Test: log-likelihood functions
# ============================================================================ #

test_that("brs_loglik returns finite scalar", {
  sim <- sim_fixed()
  ll <- betaregscale:::brs_loglik(
    param = c(0.2, -0.5, 0.3, 1 / 5),
    formula = y ~ x1 + x2, data = sim
  )
  expect_true(is.numeric(ll))
  expect_equal(length(ll), 1L)
  expect_true(is.finite(ll))
  expect_true(ll < 0)
})

test_that("internal variable loglik wrapper returns finite scalar", {
  sim <- sim_variable()
  ll <- betaregscale:::brs_loglik_var(
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

test_that("brs_fit_fixed converges and returns correct class", {
  sim <- sim_fixed()
  fit <- brs_fit_fixed(
    formula = y ~ x1 + x2, data = sim,
    link = "logit", link_phi = "logit"
  )
  expect_s3_class(fit, "brs")
  expect_equal(fit$convergence, 0L)
  expect_true(is.finite(fit$value))
  expect_equal(length(fit$par), 4L)
  expect_equal(fit$nobs, nrow(sim))
})

test_that("brs_fit_var converges for variable dispersion", {
  sim <- sim_variable()
  fit <- brs_fit_var(
    formula = y ~ x1 + x2 | z1, data = sim,
    link = "logit", link_phi = "logit"
  )
  expect_s3_class(fit, "brs")
  expect_equal(fit$convergence, 0L)
  expect_equal(length(fit$par), 5L)
})

test_that("brs() dispatches to fixed model", {
  sim <- sim_fixed()
  fit <- brs(y ~ x1 + x2, data = sim)
  expect_s3_class(fit, "brs")
  expect_equal(length(fit$par), 4L)
})

test_that("brs() dispatches to variable model", {
  sim <- sim_variable()
  fit <- brs(y ~ x1 + x2 | z1, data = sim)
  expect_s3_class(fit, "brs")
  expect_equal(length(fit$par), 5L)
})

test_that("fitted object stores call", {
  sim <- sim_fixed()
  fit <- brs(y ~ x1 + x2, data = sim)
  expect_true(!is.null(fit$call))
})


# ============================================================================ #
# Test: parameter naming (betareg style)
# ============================================================================ #

test_that("fixed model parameters have correct names", {
  sim <- sim_fixed()
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)
  nm <- names(fit$par)
  expect_equal(nm, c("(Intercept)", "x1", "x2", "(phi)"))
  expect_equal(names(fit$coefficients$mean), c("(Intercept)", "x1", "x2"))
  expect_equal(names(fit$coefficients$precision), "(phi)")
})

test_that("variable model parameters have (phi)_ prefix", {
  sim <- sim_variable()
  fit <- brs_fit_var(y ~ x1 + x2 | z1, data = sim)
  nm <- names(fit$par)
  expect_true("(phi)_(Intercept)" %in% nm)
  expect_true("(phi)_z1" %in% nm)
})


# ============================================================================ #
# Test: S3 methods <U+2014> coef with model= argument
# ============================================================================ #

test_that("coef method works with model argument", {
  sim <- sim_fixed()
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)

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
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)

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
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_true(is.finite(as.numeric(ll)))
  expect_equal(attr(ll, "df"), 4L)
})

test_that("AIC and BIC return scalars", {
  sim <- sim_fixed()
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)
  aic <- AIC(fit)
  bic <- BIC(fit)
  expect_true(is.finite(aic))
  expect_true(is.finite(bic))
  expect_true(bic > aic) # BIC penalizes more with n > e^2
})

test_that("stats::AIC(logLik(fit)) == AIC(fit)", {
  sim <- sim_fixed()
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)
  expect_equal(stats::AIC(logLik(fit)), AIC(fit), tolerance = 1e-6)
})

test_that("stats::BIC(logLik(fit)) == BIC(fit)", {
  sim <- sim_fixed()
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)
  expect_equal(stats::BIC(logLik(fit)), BIC(fit), tolerance = 1e-6)
})


# ============================================================================ #
# Test: nobs, formula, model.matrix
# ============================================================================ #

test_that("nobs returns correct count", {
  sim <- sim_fixed()
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)
  expect_equal(nobs(fit), 200L)
})

test_that("formula returns the model formula", {
  sim <- sim_fixed()
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)
  f <- formula(fit)
  expect_true(inherits(f, "formula"))
})

test_that("model.matrix returns design matrices", {
  sim <- sim_variable()
  fit <- brs_fit_var(y ~ x1 + x2 | z1, data = sim)
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
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)
  s <- summary(fit)
  expect_s3_class(s, "summary.brs")
  expect_true("coefficients" %in% names(s))
  expect_true("mean" %in% names(s$coefficients))
  expect_true("precision" %in% names(s$coefficients))
  expect_equal(nrow(s$coefficients$mean), 3L)
  expect_equal(nrow(s$coefficients$precision), 1L)
  expect_true(!is.null(s$censoring))
})

test_that("summary uses pnorm not pt for p-values", {
  sim <- sim_fixed()
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)
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
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)
  expect_output(print(fit), "Coefficients")
})

test_that("print.summary method runs without error", {
  sim <- sim_fixed()
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)
  expect_output(print(summary(fit)), "Coefficients")
})


# ============================================================================ #
# Test: fitted values
# ============================================================================ #

test_that("fitted method returns vectors of correct length", {
  sim <- sim_fixed()
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)
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
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)

  for (rtype in c("response", "pearson", "deviance", "rqr")) {
    r <- residuals(fit, type = rtype)
    expect_equal(length(r), nrow(sim))
    expect_true(all(is.finite(r)))
  }
})

test_that("RQR residuals are approximately normal for correct model", {
  sim <- sim_fixed(n = 500, seed = 123)
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)
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
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)
  ci <- confint(fit)
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 4L)
  expect_equal(ncol(ci), 2L)
  expect_true(all(ci[, 1] < ci[, 2]))
})

test_that("confint with model argument works", {
  sim <- sim_fixed()
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)
  ci_mean <- confint(fit, model = "mean")
  ci_prec <- confint(fit, model = "precision")
  expect_equal(nrow(ci_mean), 3L)
  expect_equal(nrow(ci_prec), 1L)
})

test_that("confint with parm argument works", {
  sim <- sim_fixed()
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)
  ci <- confint(fit, parm = c("x1", "x2"))
  expect_equal(nrow(ci), 2L)
})


# ============================================================================ #
# Test: predict
# ============================================================================ #

test_that("predict with type='response' works", {
  sim <- sim_fixed()
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)
  pred <- predict(fit, type = "response")
  expect_equal(length(pred), nrow(sim))
  expect_true(all(pred > 0 & pred < 1))
})

test_that("predict with newdata works", {
  sim <- sim_fixed()
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)
  nd <- data.frame(x1 = c(-1, 0, 1), x2 = c(0, 0, 0))
  pred <- predict(fit, newdata = nd, type = "response")
  expect_equal(length(pred), 3L)
  expect_true(all(pred > 0 & pred < 1))
})

test_that("predict type='link' returns linear predictor", {
  sim <- sim_fixed()
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)
  pred <- predict(fit, type = "link")
  expect_equal(length(pred), nrow(sim))
  # link values are on unrestricted scale
})

test_that("predict type='variance' works", {
  sim <- sim_fixed()
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)
  pred <- predict(fit, type = "variance")
  expect_equal(length(pred), nrow(sim))
  expect_true(all(pred > 0))
})

test_that("predict type='quantile' works", {
  sim <- sim_fixed()
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)
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
  fit <- brs_fit_var(y ~ x1 + x2 | z1, data = sim)
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
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)
  g <- brs_gof(fit)
  expect_true(is.data.frame(g))
  expect_true(all(c("logLik", "AIC", "BIC") %in% names(g)))
})

test_that("est returns estimates with p-values using pnorm", {
  sim <- sim_fixed()
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)
  e <- brs_est(fit)
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
  sim <- brs_sim(
    formula = ~x1, data = dat,
    beta = c(0.2, 0.3), phi = 1 / 5,
    link = "logit", link_phi = "logit"
  )
  fit <- brs_fit_fixed(
    y ~ x1,
    data = sim,
    link = "probit", link_phi = "logit"
  )
  expect_s3_class(fit, "brs")
  expect_equal(fit$convergence, 0L)
})

test_that("brs_coef produces valid confidence intervals", {
  sim <- sim_fixed()
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)
  tab <- brs_coef(fit)$est
  expect_true(all(tab$ci_lower < tab$estimate))
  expect_true(all(tab$ci_upper > tab$estimate))
})


# ============================================================================ #
# Test: mixed censoring in brs_check
# ============================================================================ #

test_that("mixed censoring: all four types present", {
  # brs_prep with all censoring types
  y_mixed <- c(0, 0, 3, 5, 7, 10, 10)
  out <- brs_prep(data.frame(y = y_mixed), ncuts = 10)
  deltas <- out[, "delta"]
  expect_true(1L %in% deltas) # left-censored
  expect_true(2L %in% deltas) # right-censored
  expect_true(3L %in% deltas) # interval-censored
})


# ============================================================================ #
# Test: plot runs without error
# ============================================================================ #

test_that("plot.brs runs without error", {
  sim <- sim_fixed()
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)
  pdf(nullfile())
  on.exit(dev.off(), add = TRUE)
  expect_silent(plot(fit, which = 1:4))
})


# ============================================================================ #
# Test: censoring_summary runs without error
# ============================================================================ #

test_that("censoring_summary returns correct structure", {
  sim <- sim_fixed()
  fit <- brs_fit_fixed(y ~ x1 + x2, data = sim)
  cs <- brs_cens(fit)
  expect_true(is.data.frame(cs))
  expect_true("type" %in% names(cs))
  expect_true("count" %in% names(cs))
  expect_true("proportion" %in% names(cs))
  expect_equal(sum(cs$count), nrow(sim))
})

test_that("censoring_summary works with raw matrix", {
  y <- c(0, 3, 5, 7, 10)
  Y <- brs_check(y, ncuts = 10)
  cs <- brs_cens(Y)
  expect_true(is.data.frame(cs))
  expect_equal(sum(cs$count), 5L)
})


# ============================================================================ #
# brs_prep() tests
# ============================================================================ #

# -- Structure and attributes ------------------------------------------------ #

test_that("brs_prep output has correct structure", {
  d <- data.frame(y = c(0, 5, 10), x1 = rnorm(3))
  res <- brs_prep(d, ncuts = 10)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("left", "right", "yt", "y", "delta") %in% names(res)))
  expect_true("x1" %in% names(res))
  expect_true(isTRUE(attr(res, "is_prepared")))
  expect_equal(attr(res, "ncuts"), 10L)
  expect_null(attr(res, "type"))
  expect_equal(attr(res, "lim"), 0.5)
})

test_that("brs_prep preserves original y values", {
  d <- data.frame(y = c(10, 20, 30), x = 1:3)
  res <- brs_prep(d, ncuts = 100)
  expect_equal(res$y, c(10, 20, 30))
})

# -- Mode 1: y only (automatic classification) ------------------------------ #

test_that("brs_prep Mode 1: y only — delta inferred correctly", {
  d <- data.frame(y = c(0, 1, 5, 9, 10), x = 1:5)
  res <- brs_prep(d, ncuts = 10)
  expect_equal(res$delta, c(1L, 3L, 3L, 3L, 2L))
})

test_that("brs_prep Mode 1: y only — endpoints match brs_check", {
  # For exact/left/right, logic is same. Interval is center.
  y_vec <- c(0, 50, 100)
  cr <- brs_check(y_vec, ncuts = 100, lim = 0.5)
  d <- data.frame(y = y_vec, x = 1:3)
  res <- brs_prep(d, ncuts = 100, lim = 0.5)
  expect_equal(res$left, as.numeric(cr[, "left"]), tolerance = 1e-8)
  expect_equal(res$right, as.numeric(cr[, "right"]), tolerance = 1e-8)
  expect_equal(res$yt, as.numeric(cr[, "yt"]), tolerance = 1e-8)
  expect_equal(res$delta, as.integer(cr[, "delta"]))
})

# -- Mode 2: y + explicit delta --------------------------------------------- #

test_that("brs_prep Mode 2: y + delta explicit — exact", {
  d <- data.frame(y = 50, delta = 0, x = 1) # Forced exact
  res <- brs_prep(d, ncuts = 100)
  expect_equal(res$delta, 0L)
  expect_equal(res$left, 0.50, tolerance = 1e-4)
  expect_equal(res$right, 0.50, tolerance = 1e-4)
})

test_that("brs_prep Mode 2: y + delta explicit — left-censored", {
  d <- data.frame(y = 0, delta = 1, x = 1)
  res <- brs_prep(d, ncuts = 100)
  expect_equal(res$delta, 1L)
  expect_true(res$left < 0.01)
  expect_equal(res$right, 0.005, tolerance = 1e-4)
})

test_that("brs_prep Mode 2: y + delta explicit — right-censored", {
  d <- data.frame(y = 100, delta = 2, x = 1)
  res <- brs_prep(d, ncuts = 100)
  expect_equal(res$delta, 2L)
  expect_equal(res$left, 0.995, tolerance = 1e-4)
  expect_true(res$right > 0.99)
})

test_that("brs_prep Mode 2: y + delta explicit — interval-censored", {
  d <- data.frame(y = 50, delta = 3, x = 1)
  res <- brs_prep(d, ncuts = 100)
  expect_equal(res$delta, 3L)
  expect_equal(res$left, 0.495, tolerance = 1e-4)
  expect_equal(res$right, 0.505, tolerance = 1e-4)
  expect_equal(res$yt, 0.50, tolerance = 1e-4)
})

# -- Mode 3: left/right with NA patterns ------------------------------------ #

test_that("brs_prep Mode 3: only right given — left-censored", {
  d <- data.frame(right = 5, x = 1) # Implies left=NA -> left-censoring
  res <- brs_prep(d, ncuts = 100)
  expect_equal(res$delta, 1L)
  expect_true(res$left < 0.001)
  expect_equal(res$right, 0.05, tolerance = 1e-4)
})

test_that("brs_prep Mode 3: only left given — right-censored", {
  d <- data.frame(left = 95, x = 1) # Implies right=NA -> right-censoring
  res <- brs_prep(d, ncuts = 100)
  expect_equal(res$delta, 2L)
  expect_equal(res$left, 0.95, tolerance = 1e-4)
  expect_true(res$right > 0.99)
})

test_that("brs_prep Mode 3: left + right given — interval-censored", {
  d <- data.frame(left = 40, right = 60, x = 1)
  res <- brs_prep(d, ncuts = 100)
  expect_equal(res$delta, 3L)
  expect_equal(res$left, 0.40, tolerance = 1e-4)
  expect_equal(res$right, 0.60, tolerance = 1e-4)
  expect_equal(res$yt, 0.50, tolerance = 1e-4)
})

test_that("brs_prep Mode 3: left=NA, right=NA, y=value -> exact (delta=0)", {
  # When analyst provides left/right columns but both are NA for a row,

  # while y has a value, the observation is exact (no censoring)
  d <- data.frame(left = NA_real_, right = NA_real_, y = 50)
  res <- brs_prep(d, ncuts = 100)
  expect_equal(res$delta, 0L)
  expect_equal(res$yt, 0.50, tolerance = 1e-4)
})

test_that("brs_prep: continuous y in (0,1) → exact", {
  d <- data.frame(y = c(0.2, 0.5, 0.8))
  res <- brs_prep(d, ncuts = 100)
  expect_equal(res$delta, c(0L, 0L, 0L))
  expect_equal(res$yt, c(0.2, 0.5, 0.8), tolerance = 1e-4)
})

# -- Mode 4: y + left + right (analyst-supplied intervals) ------------------- #

test_that("brs_prep Mode 4: y + left + right — uses analyst endpoints", {
  d <- data.frame(y = 50, left = 48, right = 52, x = 1)
  res <- brs_prep(d, ncuts = 100)
  expect_equal(res$delta, 3L)
  expect_equal(res$left, 0.48, tolerance = 1e-4)
  expect_equal(res$right, 0.52, tolerance = 1e-4)
  expect_equal(res$yt, 0.50, tolerance = 1e-4)
})

# -- Mixed censoring in one dataset ----------------------------------------- #

test_that("brs_prep handles mixed censoring types", {
  d <- data.frame(
    y = c(0, 100, 50, NA),
    left = c(NA, NA, NA, 30),
    right = c(NA, NA, NA, 70),
    x = 1:4
  )
  res <- brs_prep(d, ncuts = 100)
  # row 4: left=NA, right=NA, y=50 -> exact (delta=0)
  expect_equal(res$delta, c(1L, 2L, 0L, 3L))
})

# -- Input validation errors ------------------------------------------------- #

test_that("brs_prep errors on non-data.frame", {
  expect_error(brs_prep(list(y = 1:5)), "data.frame")
})

test_that("brs_prep errors when no relevant columns", {
  expect_error(brs_prep(data.frame(x1 = 1:5)), "must be present")
})

test_that("brs_prep errors on invalid delta values", {
  d <- data.frame(y = 50, delta = 99)
  expect_error(brs_prep(d, ncuts = 100), "\\{0, 1, 2, 3\\}")
})

test_that("brs_prep errors when left > right", {
  d <- data.frame(left = 60, right = 40)
  expect_error(brs_prep(d, ncuts = 100), "left > right")
})

test_that("brs_prep errors when all columns are NA", {
  d <- data.frame(y = NA_real_, left = NA_real_, right = NA_real_)
  expect_error(brs_prep(d, ncuts = 100), "all relevant columns are NA")
})

test_that("brs_prep errors when ncuts < max value", {
  d <- data.frame(y = c(50, 150))
  expect_error(brs_prep(d, ncuts = 100), "ncuts")
})

test_that("brs_prep errors on non-numeric y", {
  d <- data.frame(y = c("a", "b"))
  expect_error(brs_prep(d, ncuts = 100), "numeric")
})

test_that("brs_prep errors on negative y", {
  d <- data.frame(y = c(-1, 5))
  expect_error(brs_prep(d, ncuts = 100), "non-negative")
})

# -- Consistency warnings ---------------------------------------------------- #

capture_warnings <- function(expr) {
  warnings <- character(0)
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = warnings)
}

test_that("brs_prep warns: delta=1 but y != 0", {
  d <- data.frame(y = 5, delta = 1L)
  out <- capture_warnings(brs_prep(d, ncuts = 100))
  expect_equal(length(out$warnings), 1L)
  expect_true(grepl("delta = 1 \\(left-censored\\) but y != 0", out$warnings[1]))
})

test_that("brs_prep warns: delta=2 but y != ncuts", {
  d <- data.frame(y = 50, delta = 2L)
  out <- capture_warnings(brs_prep(d, ncuts = 100))
  expect_equal(length(out$warnings), 1L)
  expect_true(grepl("delta = 2 \\(right-censored\\) but y != 100", out$warnings[1]))
})

test_that("brs_prep warns: delta=3 but y at boundary", {
  d <- data.frame(y = 0, delta = 3L)
  out <- capture_warnings(brs_prep(d, ncuts = 100))
  expect_equal(length(out$warnings), 1L)
  expect_true(grepl("delta = 3 \\(interval-censored\\) but y is at a boundary", out$warnings[1]))
})

# -- censoring_summary with brs_prep output -------------------------------- #

test_that("censoring_summary works with brs_prep output", {
  d <- data.frame(y = c(0, 3, 5, 7, 10), x1 = rnorm(5))
  res <- brs_prep(d, ncuts = 10)
  cs <- brs_cens(res)
  expect_true(is.data.frame(cs))
  expect_equal(sum(cs$count), 5L)
})

# -- End-to-end: brs_prep -> betaregscale ---------------------------------- #

test_that("brs_prep -> betaregscale produces identical fit to raw data", {
  sim <- sim_fixed(n = 200, seed = 42)

  # Fit 1: standard workflow
  fit1 <- brs(y ~ x1 + x2, data = sim, link = "logit", repar = 2L)

  # Fit 2: brs_prep workflow <U+2014> extract only y + covariates
  sim_raw <- sim[, c("y", "x1", "x2")]
  prep <- brs_prep(sim_raw, ncuts = 100)
  fit2 <- brs(y ~ x1 + x2, data = prep, link = "logit", repar = 2L)

  expect_equal(coef(fit1), coef(fit2), tolerance = 1e-6)
  expect_equal(logLik(fit1), logLik(fit2), tolerance = 1e-6)
})

test_that("brs_prep -> brs_fit_var converges", {
  sim <- sim_variable(n = 200, seed = 42)
  sim_raw <- sim[, c("y", "x1", "x2", "z1")]
  prep <- brs_prep(sim_raw, ncuts = 100)
  fit <- brs(y ~ x1 + x2 | z1, data = prep, link = "logit", repar = 2L)
  expect_s3_class(fit, "brs")
  expect_equal(fit$convergence, 0L)
})


# ============================================================================ #
# Test: delta parameter in simulation functions
# ============================================================================ #

test_that("simulate delta=0 produces all exact observations", {
  set.seed(42)
  n <- 100
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  sim <- brs_sim(
    formula = ~ x1 + x2, data = dat,
    beta = c(0.2, -0.5, 0.3), phi = 1 / 5,
    delta = 0
  )
  expect_true(all(sim$delta == 0))
  # Exact: left == right == yt
  expect_equal(sim$left, sim$right)
  expect_equal(sim$left, sim$yt)
  # y should be in (0, 1) for exact
  expect_true(all(sim$y > 0 & sim$y < 1))
})

test_that("simulate delta=1 produces all left-censored observations", {
  set.seed(42)
  n <- 100
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  sim <- brs_sim(
    formula = ~ x1 + x2, data = dat,
    beta = c(0.2, -0.5, 0.3), phi = 1 / 5,
    delta = 1
  )
  expect_true(all(sim$delta == 1))
  # y values preserve simulated variation (not all forced to 0)
  expect_true(length(unique(sim$y)) > 1)
  # Left-censored: left = eps, right varies per observation
  expect_true(all(sim$left == 1e-5))
  expect_true(length(unique(sim$right)) > 1)
})

test_that("simulate delta=2 produces all right-censored observations", {
  set.seed(42)
  n <- 100
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  sim <- brs_sim(
    formula = ~ x1 + x2, data = dat,
    beta = c(0.2, -0.5, 0.3), phi = 1 / 5,
    delta = 2, ncuts = 100L
  )
  expect_true(all(sim$delta == 2))
  # y values preserve simulated variation (not all forced to 100)
  expect_true(length(unique(sim$y)) > 1)
  # Right-censored: right = 1-eps, left varies per observation
  expect_true(all(sim$right == 1 - 1e-5))
  expect_true(length(unique(sim$left)) > 1)
})

test_that("simulate delta=3 produces all interval-censored observations", {
  set.seed(42)
  n <- 100
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  sim <- brs_sim(
    formula = ~ x1 + x2, data = dat,
    beta = c(0.2, -0.5, 0.3), phi = 1 / 5,
    delta = 3, ncuts = 100L
  )
  expect_true(all(sim$delta == 3))
  # No boundary values
  expect_true(all(sim$y > 0 & sim$y < 100))
})

test_that("simulate delta=NULL preserves default mixed behavior", {
  set.seed(42)
  n <- 200
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  sim <- brs_sim(
    formula = ~ x1 + x2, data = dat,
    beta = c(0.2, -0.5, 0.3), phi = 1 / 5
  )
  # Should have a mix of censoring types (at least interval)
  expect_true(3L %in% sim$delta)
})

test_that("simulate errors on invalid delta", {
  set.seed(42)
  n <- 10
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  expect_error(
    brs_sim(
      formula = ~ x1 + x2, data = dat,
      beta = c(0.2, -0.5, 0.3), phi = 1 / 5,
      delta = 5
    ),
    "\\{0, 1, 2, 3\\}"
  )
  expect_error(
    brs_sim(
      formula = ~ x1 + x2, data = dat,
      beta = c(0.2, -0.5, 0.3), phi = 1 / 5,
      delta = c(0, 1)
    ),
    "single integer"
  )
})

test_that("simulate_z delta=3 produces all interval-censored", {
  set.seed(42)
  n <- 100
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n), z1 = runif(n))
  sim <- brs_sim(
    formula = ~ x1 + x2 | z1, data = dat,
    beta = c(0.2, -0.5, 0.3), zeta = c(0.5, -0.5),
    delta = 3
  )
  expect_true(all(sim$delta == 3))
  expect_true(all(sim$y > 0 & sim$y < 100))
})

test_that("simulate_z delta=0 produces all exact observations", {
  set.seed(42)
  n <- 100
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n), z1 = runif(n))
  sim <- brs_sim(
    formula = ~ x1 + x2 | z1, data = dat,
    beta = c(0.2, -0.5, 0.3), zeta = c(0.5, -0.5),
    delta = 0
  )
  expect_true(all(sim$delta == 0))
  expect_equal(sim$left, sim$right)
})

# -- End-to-end: simulate with forced delta -> fit ----------------------------- #

test_that("simulate delta=3 -> betaregscale converges", {
  set.seed(42)
  n <- 200
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  sim <- brs_sim(
    formula = ~ x1 + x2, data = dat,
    beta = c(0.2, -0.5, 0.3), phi = 1 / 5,
    delta = 3
  )
  fit <- brs(y ~ x1 + x2, data = sim)
  expect_equal(fit$convergence, 0L)
})

test_that("simulate delta=0 -> betaregscale converges", {
  set.seed(42)
  n <- 200
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  sim <- brs_sim(
    formula = ~ x1 + x2, data = dat,
    beta = c(0.2, -0.5, 0.3), phi = 1 / 5,
    delta = 0
  )
  fit <- brs(y ~ x1 + x2, data = sim)
  expect_equal(fit$convergence, 0L)
})

test_that("simulate delta=1 -> betaregscale converges", {
  set.seed(42)
  n <- 200
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  sim <- brs_sim(
    formula = ~ x1 + x2, data = dat,
    beta = c(0.2, -0.5, 0.3), phi = 1 / 5,
    delta = 1
  )
  fit <- brs(y ~ x1 + x2, data = sim)
  expect_equal(fit$convergence, 0L)
})

test_that("simulate delta=2 -> betaregscale converges", {
  set.seed(42)
  n <- 200
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  sim <- brs_sim(
    formula = ~ x1 + x2, data = dat,
    beta = c(0.2, -0.5, 0.3), phi = 1 / 5,
    delta = 2
  )
  fit <- brs(y ~ x1 + x2, data = sim)
  expect_equal(fit$convergence, 0L)
})

# ============================================================================ #
# Test: brs_check with user-supplied delta
# ============================================================================ #

test_that("brs_check with delta vector overrides boundary rules", {
  y <- c(0, 50, 100)
  d <- c(0L, 0L, 0L) # force all to exact
  res <- brs_check(y, ncuts = 100, delta = d)
  expect_true(all(res[, "delta"] == 0))
})

test_that("brs_check with delta vector respects user classification", {
  y <- c(50, 50)
  d <- c(1L, 2L)
  res <- brs_check(y, ncuts = 100, delta = d)
  expect_equal(as.integer(res[, "delta"]), c(1L, 2L))
})

test_that("brs_check forced delta=1 on non-boundary uses y-based upper bound", {
  y <- c(30, 60)
  d <- c(1L, 1L)
  res <- brs_check(y, ncuts = 100, delta = d)
  # Left-censored: left = eps, right = (y + lim) / ncuts
  expect_equal(unname(res[, "left"]), c(1e-5, 1e-5))
  expect_equal(unname(res[, "right"]), c((30 + 0.5) / 100, (60 + 0.5) / 100),
    tolerance = 1e-6
  )
})

test_that("brs_check forced delta=2 on non-boundary uses y-based lower bound", {
  y <- c(30, 60)
  d <- c(2L, 2L)
  res <- brs_check(y, ncuts = 100, delta = d)
  # Right-censored: left = (y - lim) / ncuts, right = 1 - eps
  expect_equal(unname(res[, "left"]), c((30 - 0.5) / 100, (60 - 0.5) / 100),
    tolerance = 1e-6
  )
  expect_equal(unname(res[, "right"]), c(1 - 1e-5, 1 - 1e-5))
})

test_that("brs_check boundary delta=1 (y=0) retains original endpoint formula", {
  y <- c(0)
  d <- c(1L)
  res <- brs_check(y, ncuts = 100, delta = d)
  expect_equal(unname(res[, "left"]), 1e-5)
  expect_equal(unname(res[, "right"]), 0.5 / 100, tolerance = 1e-6)
})

test_that("brs_check boundary delta=2 (y=ncuts) retains original endpoint formula", {
  y <- c(100)
  d <- c(2L)
  res <- brs_check(y, ncuts = 100, delta = d)
  expect_equal(unname(res[, "left"]), (100 - 0.5) / 100, tolerance = 1e-6)
  expect_equal(unname(res[, "right"]), 1 - 1e-5)
})

test_that("brs_check delta must match y length", {
  y <- c(1, 2, 3)
  expect_error(brs_check(y, ncuts = 100, delta = c(0L, 1L)), "same length")
})

test_that("brs_check delta must be valid values", {
  y <- c(1, 2, 3)
  expect_error(brs_check(y, ncuts = 100, delta = c(0L, 1L, 5L)), "\\{0, 1, 2, 3\\}")
})
