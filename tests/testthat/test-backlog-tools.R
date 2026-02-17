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

test_that("brs_marginaleffects returns class and stores draws when requested", {
  sim <- .sim_for_tools(seed = 1919)
  fit <- brs(y ~ x1 + x2 | z1, data = sim, repar = 2)

  me <- brs_marginaleffects(
    fit,
    model = "mean",
    type = "response",
    interval = TRUE,
    n_sim = 80,
    keep_draws = TRUE
  )

  expect_s3_class(me, "brs_marginaleffects")
  expect_true(is.matrix(attr(me, "ame_draws")))
  expect_equal(ncol(attr(me, "ame_draws")), nrow(me))
  expect_equal(attr(me, "n_sim"), 80L)
})

test_that("autoplot.brs_marginaleffects supports forest/magnitude/dist", {
  skip_if_not_installed("ggplot2")
  sim <- .sim_for_tools(seed = 2101)
  fit <- brs(y ~ x1 + x2 | z1, data = sim, repar = 2)

  me <- brs_marginaleffects(
    fit,
    model = "mean",
    type = "response",
    interval = TRUE,
    n_sim = 80,
    keep_draws = TRUE
  )

  p1 <- autoplot.brs_marginaleffects(me, type = "forest")
  p2 <- autoplot.brs_marginaleffects(me, type = "magnitude")
  p3 <- autoplot.brs_marginaleffects(me, type = "dist", variable = me$variable[1L])

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
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

test_that("ggplot2::autoplot dispatches for brs and brs_marginaleffects", {
  skip_if_not_installed("ggplot2")
  sim <- .sim_for_tools(seed = 3030)
  fit <- brs(y ~ x1 + x2 | z1, data = sim, repar = 2)

  p_brs <- ggplot2::autoplot(fit, type = "calibration")
  expect_s3_class(p_brs, "ggplot")

  me <- brs_marginaleffects(
    fit,
    model = "mean",
    type = "response",
    interval = TRUE,
    n_sim = 80,
    keep_draws = TRUE
  )
  p_me <- ggplot2::autoplot(me, type = "forest")
  expect_s3_class(p_me, "ggplot")
})

test_that("brs_predict_scoreprob returns coherent probability matrices", {
  sim <- .sim_for_tools(seed = 404)
  fit <- brs(y ~ x1 + x2 | z1, data = sim, repar = 2)

  P <- brs_predict_scoreprob(fit)
  expect_true(is.matrix(P))
  expect_equal(nrow(P), nrow(sim))
  expect_equal(ncol(P), fit$ncuts + 1L)
  expect_true(all(rowSums(P) > 0.98 & rowSums(P) < 1.02))

  L <- brs_predict_scoreprob(fit, scores = 0:10, format = "long")
  expect_true(is.data.frame(L))
  expect_true(all(c("id", "score", "prob") %in% names(L)))
})

test_that("brs_cv returns fold-level predictive metrics", {
  sim <- .sim_for_tools(seed = 505)

  cv <- brs_cv(
    y ~ x1 + x2 | z1,
    data = sim,
    k = 3,
    repeats = 1,
    repar = 2,
    seed = 505
  )

  expect_s3_class(cv, "brs_cv")
  expect_true(is.data.frame(cv))
  expect_equal(nrow(cv), 3L)
  expect_true(all(c("repeat", "fold", "log_score", "rmse_yt", "mae_yt") %in% names(cv)))
  expect_true(all(cv$fold %in% 1:3))
})
