# ============================================================================ #
# Tests for brsmm mixed-effects models
# ============================================================================ #

.sim_brsmm <- function(seed = 2026L, g = 12L, ni = 8L) {
  set.seed(seed)
  id <- factor(rep(seq_len(g), each = ni))
  n <- length(id)
  x1 <- rnorm(n)
  b <- rnorm(g, sd = 0.45)

  eta_mu <- 0.2 + 0.6 * x1 + b[as.integer(id)]
  mu <- plogis(eta_mu)
  phi <- plogis(-0.3 + 0.2 * x1)
  shp <- brs_repar(mu = mu, phi = phi, repar = 2L)

  y <- round(stats::rbeta(n, shp$shape1, shp$shape2) * 100)
  data.frame(y = y, x1 = x1, id = id)
}

test_that("brsmm fits random-intercept models with Laplace", {
  d <- .sim_brsmm()
  fit <- brsmm(
    y ~ x1,
    random = ~ 1 | id,
    data = d,
    repar = 2L,
    int_method = "laplace",
    method = "BFGS",
    control = list(maxit = 500L)
  )

  expect_s3_class(fit, "brsmm")
  expect_true(is.finite(fit$value))
  expect_equal(fit$ngroups, nlevels(d$id))
  expect_true(all(c("mean", "precision", "random") %in% names(fit$coefficients)))
  expect_true(length(fit$random$mode_b) == nlevels(d$id))
  expect_true(is.finite(fit$random$sigma_b))
})

test_that("brsmm supports two-part formulas and prediction", {
  d <- .sim_brsmm(seed = 2027L)
  fit <- brsmm(
    y ~ x1 | x1,
    random = ~ 1 | id,
    data = d,
    repar = 2L,
    int_method = "laplace",
    control = list(maxit = 500L)
  )

  p1 <- predict(fit, type = "response")
  p2 <- predict(fit, type = "precision")
  p3 <- predict(fit, type = "variance")

  expect_equal(length(p1), nrow(d))
  expect_equal(length(p2), nrow(d))
  expect_equal(length(p3), nrow(d))
  expect_true(all(is.finite(p1)))
  expect_true(all(is.finite(p2)))
  expect_true(all(is.finite(p3)))
})

test_that("brsmm S3 methods return coherent outputs", {
  d <- .sim_brsmm(seed = 2028L, g = 10L, ni = 7L)
  fit <- brsmm(
    y ~ x1 | x1,
    random = ~ 1 | id,
    data = d,
    repar = 2L,
    int_method = "laplace",
    method = "BFGS",
    control = list(maxit = 700L)
  )

  # Coefficients by component
  cf_full <- coef(fit)
  cf_mean <- coef(fit, model = "mean")
  cf_prec <- coef(fit, model = "precision")
  cf_rand <- coef(fit, model = "random")
  expect_true(length(cf_full) == (length(cf_mean) + length(cf_prec) + length(cf_rand)))
  expect_true(all(is.finite(cf_full)))

  # Variance-covariance blocks
  v_full <- vcov(fit, model = "full")
  v_mean <- vcov(fit, model = "mean")
  v_prec <- vcov(fit, model = "precision")
  v_rand <- vcov(fit, model = "random")
  expect_equal(dim(v_full), c(length(cf_full), length(cf_full)))
  expect_equal(dim(v_mean), c(length(cf_mean), length(cf_mean)))
  expect_equal(dim(v_prec), c(length(cf_prec), length(cf_prec)))
  expect_equal(dim(v_rand), c(length(cf_rand), length(cf_rand)))

  # Likelihood and information criteria
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_equal(attr(ll, "nobs"), nrow(d))
  expect_true(is.finite(as.numeric(ll)))
  expect_true(is.finite(AIC(fit)))
  expect_true(is.finite(BIC(fit)))
  expect_equal(nobs(fit), nrow(d))

  # Fitted/predicted values
  mu_fit <- fitted(fit, type = "mu")
  phi_fit <- fitted(fit, type = "phi")
  expect_equal(length(mu_fit), nrow(d))
  expect_equal(length(phi_fit), nrow(d))
  expect_true(all(is.finite(mu_fit)))
  expect_true(all(is.finite(phi_fit)))

  p_resp <- predict(fit, type = "response")
  p_link <- predict(fit, type = "link")
  p_prec <- predict(fit, type = "precision")
  p_var <- predict(fit, type = "variance")
  expect_equal(length(p_resp), nrow(d))
  expect_equal(length(p_link), nrow(d))
  expect_equal(length(p_prec), nrow(d))
  expect_equal(length(p_var), nrow(d))
  expect_true(all(is.finite(p_resp)))
  expect_true(all(is.finite(p_link)))
  expect_true(all(is.finite(p_prec)))
  expect_true(all(is.finite(p_var)))

  # Residuals
  r_resp <- residuals(fit, type = "response")
  r_pear <- residuals(fit, type = "pearson")
  expect_equal(length(r_resp), nrow(d))
  expect_equal(length(r_pear), nrow(d))
  expect_true(all(is.finite(r_resp)))
  expect_true(all(is.finite(r_pear)))

  # Summary and printing
  sm <- summary(fit)
  expect_s3_class(sm, "summary.brsmm")
  expect_equal(nrow(sm$coefficients), length(cf_full))
  expect_true(all(c("Estimate", "Std. Error", "z value", "Pr(>|z|)") %in% colnames(sm$coefficients)))
  expect_output(print(fit), "Mixed beta interval model")
  expect_output(print(sm), "Mixed beta interval model")

  # Random effects
  re <- ranef.brsmm(fit)
  expect_true(length(re) == nlevels(d$id))
  expect_true(all(is.finite(re)))
})

test_that("brsmm predict() with newdata handles known and unseen groups", {
  d <- .sim_brsmm(seed = 2029L, g = 9L, ni = 8L)
  fit <- brsmm(
    y ~ x1 | x1,
    random = ~ 1 | id,
    data = d,
    repar = 2L,
    int_method = "laplace",
    control = list(maxit = 700L)
  )

  new_known <- d[1:6, c("x1", "id"), drop = FALSE]
  pred_known <- predict(fit, newdata = new_known, type = "response")
  expect_equal(length(pred_known), nrow(new_known))
  expect_true(all(is.finite(pred_known)))

  new_unseen <- new_known
  new_unseen$id <- factor(rep("new_group", nrow(new_unseen)))
  pred_unseen <- predict(fit, newdata = new_unseen, type = "response")
  expect_equal(length(pred_unseen), nrow(new_unseen))
  expect_true(all(is.finite(pred_unseen)))

  # Without group column in newdata, method should default to random effect = 0.
  new_nogroup <- new_known["x1"]
  pred_nogroup <- predict(fit, newdata = new_nogroup, type = "response")
  expect_equal(length(pred_nogroup), nrow(new_nogroup))
  expect_true(all(is.finite(pred_nogroup)))
})

.sim_brsmm_recovery <- function(seed = 2030L, g = 40L, ni = 10L,
                                beta = c(0.25, 0.70),
                                gamma = -0.20,
                                sigma_b = 0.55) {
  set.seed(seed)
  id <- factor(rep(seq_len(g), each = ni))
  n <- length(id)
  x1 <- rnorm(n)
  b_true <- rnorm(g, mean = 0, sd = sigma_b)

  eta_mu <- beta[1L] + beta[2L] * x1 + b_true[as.integer(id)]
  eta_phi <- rep(gamma, n)
  mu <- plogis(eta_mu)
  phi <- plogis(eta_phi)
  shp <- brs_repar(mu = mu, phi = phi, repar = 2L)

  y <- round(stats::rbeta(n, shp$shape1, shp$shape2) * 100)
  list(
    data = data.frame(y = y, x1 = x1, id = id),
    beta_true = beta,
    gamma_true = gamma,
    sigma_b_true = sigma_b,
    b_true = stats::setNames(b_true, levels(id))
  )
}

test_that("brsmm recovers fixed and random-effect scale parameters", {
  sim <- .sim_brsmm_recovery()
  d <- sim$data
  fit <- brsmm(
    y ~ x1,
    random = ~ 1 | id,
    data = d,
    repar = 2L,
    int_method = "laplace",
    method = "BFGS",
    control = list(maxit = 900L)
  )

  b_hat <- coef(fit, model = "mean")
  sigma_hat <- exp(unname(coef(fit, model = "random")))
  re_hat <- ranef.brsmm(fit)
  re_true <- sim$b_true[names(re_hat)]

  expect_lt(abs(unname(b_hat[1L]) - sim$beta_true[1L]), 0.35)
  expect_lt(abs(unname(b_hat[2L]) - sim$beta_true[2L]), 0.35)
  expect_lt(abs(log(sigma_hat) - log(sim$sigma_b_true)), 0.50)

  # Posterior modes are shrunk but should track true random effects.
  expect_gt(stats::cor(as.numeric(re_hat), as.numeric(re_true)), 0.45)
})

test_that("brsmm plot and autoplot methods run and return expected objects", {
  d <- .sim_brsmm(seed = 2031L, g = 8L, ni = 7L)
  fit <- brsmm(
    y ~ x1 | x1,
    random = ~ 1 | id,
    data = d,
    repar = 2L,
    int_method = "laplace",
    control = list(maxit = 700L)
  )

  # Base and gg backends should run and invisibly return the fitted object.
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)

  p_base <- plot(fit, which = 1:2, type = "pearson")
  expect_s3_class(p_base, "brsmm")

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p_gg <- plot(fit, which = 1:2, gg = TRUE)
    expect_s3_class(p_gg, "brsmm")

    ap1 <- autoplot.brsmm(fit, type = "calibration")
    ap2 <- autoplot.brsmm(fit, type = "score_dist")
    ap3 <- autoplot.brsmm(fit, type = "ranef_qq")
    ap4 <- autoplot.brsmm(fit, type = "residuals_by_group")

    expect_s3_class(ap1, "ggplot")
    expect_s3_class(ap2, "ggplot")
    expect_s3_class(ap3, "ggplot")
    expect_s3_class(ap4, "ggplot")
  }
})
