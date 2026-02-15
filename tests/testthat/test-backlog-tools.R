# ============================================================================ #
# Tests for new analyst-oriented utilities
# ============================================================================ #

.sim_for_tools <- function(n = 120L, seed = 123) {
  set.seed(seed)
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), z1 = rnorm(n))
  brs_sim(
    formula = ~ x1 + x2 | z1,
    data = d,
    beta = c(0.2, -0.4, 0.2),
    zeta = c(0.1, -0.2),
    ncuts = 100,
    repar = 2
  )
}

test_that("brs_table compares multiple fitted models", {
  sim <- .sim_for_tools()
  m1 <- brs(y ~ x1 + x2, data = sim, repar = 2)
  m2 <- brs(y ~ x1 + x2 | z1, data = sim, repar = 2)

  tab <- brs_table(fixed = m1, variable = m2, sort_by = "AIC")
  expect_true(is.data.frame(tab))
  expect_equal(nrow(tab), 2L)
  expect_true(all(c("model", "AIC", "BIC", "logLik", "pseudo_r2") %in% names(tab)))
  expect_true(all(c("exact", "left", "right", "interval") %in% names(tab)))
})

test_that("brs_table accepts list input and no censoring columns", {
  sim <- .sim_for_tools(seed = 777)
  m1 <- brs(y ~ x1 + x2, data = sim, repar = 2)
  m2 <- brs(y ~ x1 + x2 | z1, data = sim, repar = 2)

  tab <- brs_table(models = list(m1 = m1, m2 = m2), include_censoring = FALSE)
  expect_true(is.data.frame(tab))
  expect_equal(nrow(tab), 2L)
  expect_false(any(c("exact", "left", "right", "interval") %in% names(tab)))
})

test_that("brs_marginaleffects returns AME table for mean model", {
  sim <- .sim_for_tools(seed = 303)
  fit <- brs(y ~ x1 + x2 | z1, data = sim, repar = 2)

  me <- brs_marginaleffects(
    fit,
    model = "mean",
    type = "response",
    interval = TRUE,
    n_sim = 120
  )
  expect_true(is.data.frame(me))
  expect_true(all(c("variable", "ame", "std.error", "ci.lower", "ci.upper") %in% names(me)))
  expect_true(all(c("x1", "x2") %in% me$variable))
})

test_that("brs_marginaleffects supports precision effects", {
  sim <- .sim_for_tools(seed = 909)
  fit <- brs(y ~ x1 + x2 | z1, data = sim, repar = 2)

  me <- brs_marginaleffects(
    fit,
    model = "precision",
    type = "link",
    interval = FALSE
  )
  expect_true(is.data.frame(me))
  expect_true("z1" %in% me$variable)
  expect_true(all(is.na(me$std.error)))
})

test_that("autoplot.brs supports calibration and score_dist", {
  skip_if_not_installed("ggplot2")
  sim <- .sim_for_tools(seed = 1201)
  fit <- brs(y ~ x1 + x2 | z1, data = sim, repar = 2)

  p1 <- autoplot.brs(fit, type = "calibration")
  p2 <- autoplot.brs(fit, type = "score_dist")

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("autoplot.brs supports cdf and residuals_by_delta", {
  skip_if_not_installed("ggplot2")
  sim <- .sim_for_tools(seed = 2002)
  fit <- brs(y ~ x1 + x2 | z1, data = sim, repar = 2)

  p1 <- autoplot.brs(fit, type = "cdf", max_curves = 4)
  p2 <- autoplot.brs(fit, type = "residuals_by_delta", residual_type = "rqr")

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})
