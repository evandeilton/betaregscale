# Tests for brsmm integration methods
# These tests verify that the new backend methods (AGHQ, QMC) run and produce
# consistent results with the baseline Laplace approximation.

test_that("brsmm runs with different integration methods", {
  skip_on_cran()

  # 1. Simulate simple random intercept data
  set.seed(123)
  n_groups <- 10
  n_per_group <- 20
  n <- n_groups * n_per_group
  id <- factor(rep(1:n_groups, each = n_per_group))
  x <- rnorm(n)

  # Parameters
  beta <- c(0.5, 1.0)
  sigma_b <- 0.3
  precision <- 15

  b <- rnorm(n_groups, 0, sigma_b)
  eta <- beta[1] + beta[2] * x + b[id]
  mu <- plogis(eta)

  # Generate beta responses
  y_raw <- rbeta(n, mu * precision, (1 - mu) * precision)
  # Scale to (0,1) for simplicity in this test (no censoring for speed
  # unless we want to test that specific path, but backend is same)
  # Actually brsmm handles scalar y_raw by treating as exact or creating
  # trivial intervals. Let's use clean numeric for speed.
  dat <- data.frame(y = y_raw, x = x, id = id)

  # 2. Fit Laplace (Control)
  # using ncuts=10 ensures some interval censoring logic is triggered if y is continuous
  fit_lap <- brsmm(y ~ x,
    random = ~ 1 | id, data = dat,
    int_method = "laplace", ncuts = 10,
    control = list(maxit = 500)
  )

  expect_s3_class(fit_lap, "brsmm")
  expect_equal(fit_lap$int_method, "laplace")
  expect_true(fit_lap$convergence == 0)

  # 3. Fit AGHQ (k=5)
  fit_aghq <- brsmm(y ~ x,
    random = ~ 1 | id, data = dat,
    int_method = "aghq", n_points = 5, ncuts = 10,
    control = list(maxit = 500)
  )

  expect_s3_class(fit_aghq, "brsmm")
  expect_equal(fit_aghq$int_method, "aghq")
  expect_true(fit_aghq$convergence == 0)

  # Compare estimates (should be very close for simple normal Random Intercept)
  # Tolerance loosely set because AGHQ(5) is better than Laplace
  expect_equal(coef(fit_lap), coef(fit_aghq), tolerance = 0.2)

  # 4. Fit QMC (k=100)
  # QMC is stochastic/approximate, so we expect rough agreement
  fit_qmc <- brsmm(y ~ x,
    random = ~ 1 | id, data = dat,
    int_method = "qmc", qmc_points = 100, ncuts = 10,
    control = list(maxit = 500)
  )

  expect_s3_class(fit_qmc, "brsmm")
  expect_equal(fit_qmc$int_method, "qmc")
  expect_true(fit_qmc$convergence == 0)

  expect_equal(coef(fit_lap), coef(fit_qmc), tolerance = 0.3)
})

test_that("brsmm print methods show correct method", {
  skip_on_cran()

  dat <- data.frame(y = rbeta(20, 2, 2), id = factor(rep(1:2, each = 10)))

  # We just mock the object structure to avoid fitting time
  # if we can, or just fit a trivial one.
  # Let's fit trivial one.
  fit <- brsmm(y ~ 1,
    random = ~ 1 | id, data = dat, int_method = "aghq",
    control = list(maxit = 10)
  )

  out <- capture.output(print(fit))
  expect_match(paste(out, collapse = " "), "Mixed beta interval model \\(AGHQ\\)")

  out_sum <- capture.output(print(summary(fit)))
  expect_match(paste(out_sum, collapse = " "), "Mixed beta interval model \\(AGHQ\\)")
})

test_that("brsmm supports random intercept and slope with Laplace", {
  skip_on_cran()

  set.seed(321)
  n_groups <- 12
  n_per_group <- 16
  n <- n_groups * n_per_group
  id <- factor(rep(seq_len(n_groups), each = n_per_group))
  x <- rnorm(n)

  b0 <- rnorm(n_groups, sd = 0.4)
  b1 <- rnorm(n_groups, sd = 0.25)
  eta <- -0.1 + 0.7 * x + b0[id] + b1[id] * x
  mu <- plogis(eta)
  phi <- 12
  y <- rbeta(n, mu * phi, (1 - mu) * phi)
  dat <- data.frame(y = y, x = x, id = id)

  fit <- brsmm(
    y ~ x,
    random = ~ 1 + x | id,
    data = dat,
    int_method = "laplace",
    ncuts = 10,
    control = list(maxit = 800)
  )

  expect_s3_class(fit, "brsmm")
  expect_equal(fit$q_re, 2L)
  expect_equal(ncol(ranef.brsmm(fit)), 2L)
  expect_equal(length(predict(fit, type = "response")), nrow(dat))
  expect_true(is.matrix(fit$random$D))
  expect_equal(dim(fit$random$D), c(2L, 2L))
})

test_that("random-effects study and new autoplot diagnostics run", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")

  set.seed(322)
  n_groups <- 10
  n_per_group <- 14
  n <- n_groups * n_per_group
  id <- factor(rep(seq_len(n_groups), each = n_per_group))
  x <- rnorm(n)
  b0 <- rnorm(n_groups, sd = 0.35)
  b1 <- rnorm(n_groups, sd = 0.20)
  eta <- 0.1 + 0.6 * x + b0[id] + b1[id] * x
  mu <- plogis(eta)
  phi <- 11
  y <- rbeta(n, mu * phi, (1 - mu) * phi)
  dat <- data.frame(y = y, x = x, id = id)

  fit <- brsmm(
    y ~ x,
    random = ~ 1 + x | id,
    data = dat,
    int_method = "laplace",
    ncuts = 10,
    control = list(maxit = 700)
  )

  rs <- brsmm_re_study(fit)
  expect_true(inherits(rs, "brsmm_re_study"))
  expect_true(all(c("term", "sd_model", "mean_mode", "sd_mode", "shrinkage_ratio") %in% names(rs$summary)))
  expect_equal(dim(rs$D), c(2L, 2L))
  expect_equal(dim(rs$Corr), c(2L, 2L))

  p1 <- autoplot.brsmm(fit, type = "ranef_caterpillar")
  p2 <- autoplot.brsmm(fit, type = "ranef_density")
  p3 <- autoplot.brsmm(fit, type = "ranef_pairs")
  p4 <- autoplot.brsmm(fit, type = "ranef_qq")

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
})

test_that("brsmm supports random intercept and slope with aghq and qmc", {
  skip_on_cran()

  set.seed(422)
  n_groups <- 10
  n_per_group <- 15
  n <- n_groups * n_per_group
  id <- factor(rep(seq_len(n_groups), each = n_per_group))
  x <- rnorm(n)

  b0 <- rnorm(n_groups, sd = 0.3)
  b1 <- rnorm(n_groups, sd = 0.2)
  eta <- 0.2 + 0.5 * x + b0[id] + b1[id] * x
  mu <- plogis(eta)
  phi <- 15
  y <- rbeta(n, mu * phi, (1 - mu) * phi)
  dat <- data.frame(y = y, x = x, id = id)

  # Laplace base
  fit_lap <- brsmm(
    y ~ x,
    random = ~ 1 + x | id, data = dat,
    int_method = "laplace", control = list(maxit = 800)
  )

  # AGHQ multivariate (q=2)
  fit_aghq <- brsmm(
    y ~ x,
    random = ~ 1 + x | id, data = dat,
    int_method = "aghq", n_points = 3, control = list(maxit = 800)
  )

  # QMC multivariate (q=2)
  fit_qmc <- brsmm(
    y ~ x,
    random = ~ 1 + x | id, data = dat,
    int_method = "qmc", qmc_points = 100, control = list(maxit = 800)
  )

  expect_s3_class(fit_aghq, "brsmm")
  expect_equal(fit_aghq$q_re, 2L)
  expect_equal(fit_aghq$int_method, "aghq")

  expect_s3_class(fit_qmc, "brsmm")
  expect_equal(fit_qmc$q_re, 2L)
  expect_equal(fit_qmc$int_method, "qmc")

  # Coefficients should be fairly similar
  expect_equal(coef(fit_lap), coef(fit_aghq), tolerance = 0.2)
  expect_equal(coef(fit_lap), coef(fit_qmc), tolerance = 0.3)
})
