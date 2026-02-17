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
