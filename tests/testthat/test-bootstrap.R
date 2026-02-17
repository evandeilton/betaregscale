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
    c(
      "parameter", "estimate", "se_boot", "ci_lower", "ci_upper",
      "mcse_lower", "mcse_upper", "wald_lower", "wald_upper", "level"
    )
  )
  expect_equal(nrow(boot), length(fit$par))
  expect_equal(boot$parameter, names(fit$par))
  expect_equal(boot$estimate, unname(fit$par))
  expect_true(all(boot$ci_lower <= boot$estimate))
  expect_true(all(boot$ci_upper >= boot$estimate))
  expect_equal(unique(boot$level), 0.95)
  expect_equal(attr(boot, "n_success"), 25L)
  expect_equal(attr(boot, "R"), 25L)
  expect_true(attr(boot, "n_attempted") >= attr(boot, "n_success"))
  expect_equal(attr(boot, "ci_type"), "percentile")
})

test_that("brs_bootstrap print method runs without error", {
  set.seed(42)
  dat <- data.frame(x1 = rnorm(50), x2 = rnorm(50))
  sim <- brs_sim(formula = ~ x1 + x2, data = dat, beta = c(0.2, -0.5, 0.3), phi = 1 / 5)
  fit <- brs(y ~ x1 + x2, data = sim)
  boot <- brs_bootstrap(fit, R = 15L, seed = 2)
  expect_output(print(boot), "Bootstrap confidence intervals")
  expect_output(print(boot), "Successful replicates")
  expect_output(print(boot), "CI:")
  expect_output(print(boot), "Attempts:")
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
  expect_error(brs_bootstrap(fit, max_tries = 9L), ">= R")
})

test_that("brs_bootstrap rejects brsmm with clear message", {
  # Object with both classes (e.g. if brsmm ever inherited from brs)
  fake_mm <- structure(list(), class = c("brsmm", "brs"))
  expect_error(brs_bootstrap(fake_mm), "does not support 'brsmm'")
})

test_that("brs_bootstrap supports ci_type and keep_draws", {
  set.seed(123)
  dat <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  sim <- brs_sim(formula = ~ x1 + x2, data = dat, beta = c(0.1, -0.2, 0.3), phi = 1 / 5)
  fit <- brs(y ~ x1 + x2, data = sim)

  boot_basic <- brs_bootstrap(fit, R = 20L, ci_type = "basic", seed = 10)
  boot_norm <- brs_bootstrap(fit, R = 20L, ci_type = "normal", seed = 10, keep_draws = TRUE)

  expect_equal(attr(boot_basic, "ci_type"), "basic")
  expect_equal(attr(boot_norm, "ci_type"), "normal")
  expect_true(is.matrix(attr(boot_norm, "boot_draws")))
  expect_equal(ncol(attr(boot_norm, "boot_draws")), length(fit$par))
  expect_equal(nrow(attr(boot_norm, "boot_draws")), attr(boot_norm, "n_success"))
})

test_that("brs_bootstrap supports bca and Monte Carlo diagnostics", {
  set.seed(2026)
  dat <- data.frame(x1 = rnorm(40), x2 = rnorm(40))
  sim <- brs_sim(formula = ~ x1 + x2, data = dat, beta = c(0.2, -0.3, 0.4), phi = 1 / 5)
  fit <- brs(y ~ x1 + x2, data = sim)

  boot_bca <- brs_bootstrap(fit, R = 12L, ci_type = "bca", seed = 7)

  expect_equal(attr(boot_bca, "ci_type"), "bca")
  expect_true(all(c("mcse_lower", "mcse_upper", "wald_lower", "wald_upper") %in% names(boot_bca)))
  expect_true(all(is.finite(boot_bca$wald_lower)))
  expect_true(all(is.finite(boot_bca$wald_upper)))
})

test_that("autoplot.brs_bootstrap returns ggplot objects", {
  skip_if_not_installed("ggplot2")

  set.seed(321)
  dat <- data.frame(x1 = rnorm(70), x2 = rnorm(70))
  sim <- brs_sim(formula = ~ x1 + x2, data = dat, beta = c(0.2, -0.1, 0.4), phi = 1 / 5)
  fit <- brs(y ~ x1 + x2, data = sim)
  boot <- brs_bootstrap(fit, R = 20L, seed = 9, keep_draws = TRUE)

  p1 <- autoplot.brs_bootstrap(boot, type = "ci_forest")
  p2 <- autoplot.brs_bootstrap(boot, type = "dist", parameter = boot$parameter[1L])
  p3 <- autoplot.brs_bootstrap(boot, type = "qq")
  p4 <- autoplot.brs_bootstrap(boot, type = "stability")

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
})

test_that("autoplot.brs_bootstrap supports title, caption and theme", {
  skip_if_not_installed("ggplot2")

  set.seed(111)
  dat <- data.frame(x1 = rnorm(65), x2 = rnorm(65))
  sim <- brs_sim(formula = ~ x1 + x2, data = dat, beta = c(0.2, -0.1, 0.35), phi = 1 / 5)
  fit <- brs(y ~ x1 + x2, data = sim)
  boot <- brs_bootstrap(fit, R = 20L, seed = 8, keep_draws = TRUE)

  p <- autoplot.brs_bootstrap(
    boot,
    type = "ci_forest",
    title = "Titulo customizado",
    caption = list("Legenda 1", "Legenda 2", "Legenda 3", "Legenda 4"),
    theme = ggplot2::theme_bw()
  )
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.brs_bootstrap checks draws and parameter validity", {
  skip_if_not_installed("ggplot2")

  set.seed(999)
  dat <- data.frame(x1 = rnorm(55), x2 = rnorm(55))
  sim <- brs_sim(formula = ~ x1 + x2, data = dat, beta = c(0.2, -0.2, 0.3), phi = 1 / 5)
  fit <- brs(y ~ x1 + x2, data = sim)
  boot_nodraw <- brs_bootstrap(fit, R = 15L, seed = 3, keep_draws = FALSE)

  expect_error(
    autoplot.brs_bootstrap(boot_nodraw, type = "dist"),
    "keep_draws = TRUE"
  )

  boot_draw <- brs_bootstrap(fit, R = 15L, seed = 3, keep_draws = TRUE)
  expect_error(
    autoplot.brs_bootstrap(boot_draw, type = "dist", parameter = "invalid_parameter"),
    "must be one of"
  )
})
