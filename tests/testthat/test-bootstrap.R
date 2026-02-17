# ============================================================================ #
# Tests for brs_bootstrap
# ============================================================================ #

test_that("brs_bootstrap returns correct structure for fixed dispersion", {
  set.seed(42)
  n <- 60
  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  sim <- brs_sim(
    formula = ~ x1 + x2, data = dat,
    beta = c(0.2, -0.5, 0.3), phi = 1 / 5, ncuts = 50
  )
  fit <- brs(y ~ x1 + x2, data = sim)
  boot <- brs_bootstrap(fit, R = 25L, level = 0.95, seed = 1)

  expect_s3_class(boot, "brs_bootstrap")
  expect_true(is.data.frame(boot))
  expect_equal(
    colnames(boot),
    c("parameter", "estimate", "se_boot", "ci_lower", "ci_upper", "level")
  )
  expect_equal(nrow(boot), length(fit$par))
  expect_equal(boot$parameter, names(fit$par))
  expect_equal(boot$estimate, unname(fit$par))
  expect_true(all(boot$ci_lower <= boot$estimate))
  expect_true(all(boot$ci_upper >= boot$estimate))
  expect_equal(unique(boot$level), 0.95)
  expect_equal(attr(boot, "n_success"), 25L)
  expect_true(attr(boot, "R") >= 25L)
})

test_that("brs_bootstrap print method runs without error", {
  set.seed(42)
  dat <- data.frame(x1 = rnorm(50), x2 = rnorm(50))
  sim <- brs_sim(formula = ~ x1 + x2, data = dat, beta = c(0.2, -0.5, 0.3), phi = 1 / 5)
  fit <- brs(y ~ x1 + x2, data = sim)
  boot <- brs_bootstrap(fit, R = 15L, seed = 2)
  expect_output(print(boot), "Bootstrap confidence intervals")
  expect_output(print(boot), "Successful replicates")
})

test_that("brs_bootstrap errors on non-brs object", {
  expect_error(brs_bootstrap(list(a = 1)), "must be a fitted 'brs' object")
})

test_that("brs_bootstrap errors on invalid R or level", {
  set.seed(42)
  dat <- data.frame(x1 = rnorm(30), x2 = rnorm(30))
  sim <- brs_sim(formula = ~ x1 + x2, data = dat, beta = c(0.2, -0.5, 0.3), phi = 1 / 5)
  fit <- brs(y ~ x1 + x2, data = sim)
  expect_error(brs_bootstrap(fit, R = 5L), "at least 10")
  expect_error(brs_bootstrap(fit, level = 0), "in \\(0, 1\\)")
  expect_error(brs_bootstrap(fit, level = 1), "in \\(0, 1\\)")
})

test_that("brs_bootstrap rejects brsmm with clear message", {
  # Object with both classes (e.g. if brsmm ever inherited from brs)
  fake_mm <- structure(list(), class = c("brsmm", "brs"))
  expect_error(brs_bootstrap(fake_mm), "does not support 'brsmm'")
})
