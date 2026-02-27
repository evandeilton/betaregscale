# ============================================================================ #
# Tests for RE study, autoplot, and display-name improvements (today's work)
# ============================================================================ #

# ---------------------------------------------------------------------------- #
# Shared fixtures
# ---------------------------------------------------------------------------- #

.make_brs_fit <- function(seed = 42L) {
  set.seed(seed)
  dat <- data.frame(
    y = c(
      0, 5, 20, 50, 75, 90, 100, 30, 60, 45,
      10, 40, 55, 70, 85, 25, 35, 65, 80, 15
    ),
    x1 = rep(c(1, 2), 10)
  )
  prep <- brs_prep(dat, ncuts = 100)
  brs(y ~ x1, data = prep)
}

.make_brsmm_fit <- function(seed = 2040L, g = 10L, ni = 6L) {
  set.seed(seed)
  id <- factor(rep(seq_len(g), each = ni))
  n <- length(id)
  x1 <- rnorm(n)
  b <- rnorm(g, sd = 0.5)
  eta <- 0.3 + 0.5 * x1 + b[as.integer(id)]
  mu <- plogis(eta)
  phi <- plogis(rep(-0.2, n))
  shp <- brs_repar(mu = mu, phi = phi, repar = 2L)
  y <- round(stats::rbeta(n, shp$shape1, shp$shape2) * 100)
  d <- data.frame(y = y, x1 = x1, id = id)
  prep <- brs_prep(d, ncuts = 100)
  brsmm(y ~ x1,
    random = ~ 1 | id, data = prep,
    repar = 2L, int_method = "laplace",
    method = "BFGS", control = list(maxit = 800L)
  )
}


# ---------------------------------------------------------------------------- #
# A: Display-name helpers (.pretty_phi_names / .pretty_re_names)
# ---------------------------------------------------------------------------- #

test_that(".pretty_phi_names strips (phi)_ prefix correctly", {
  nms <- c("(phi)_(Intercept)", "(phi)_x1", "other_name")
  out <- betaregscale:::.pretty_phi_names(nms)
  expect_equal(out, c("(Intercept)", "x1", "other_name"))
})

test_that(".pretty_re_names converts internal Cholesky names", {
  nms <- c(
    "(re_chol_logsd)_(Intercept)|id",
    "(re_chol)_(Intercept):x1|id",
    "plain_name"
  )
  out <- betaregscale:::.pretty_re_names(nms)
  expect_equal(out, c("logSD.(Intercept)|id", "cov.(Intercept):x1|id", "plain_name"))
})

test_that("print.summary.brsmm uses cleaned phi labels in output", {
  fit <- .make_brsmm_fit()
  sm <- summary(fit)
  out <- capture.output(print(sm))
  expect_false(any(grepl("\\(phi\\)_", out)))
})


# ---------------------------------------------------------------------------- #
# B: brsmm_re_study — ICC and VarCorr
# ---------------------------------------------------------------------------- #

test_that("brsmm_re_study returns icc in [0, 1]", {
  fit <- .make_brsmm_fit()
  rs <- brsmm_re_study(fit)
  expect_s3_class(rs, "brsmm_re_study")
  expect_true("icc" %in% names(rs))
  expect_gte(rs$icc, 0)
  expect_lte(rs$icc, 1)
})

test_that("print.brsmm_re_study shows VarCorr block and ICC", {
  fit <- .make_brsmm_fit()
  rs <- brsmm_re_study(fit)
  out <- capture.output(print(rs))
  expect_true(any(grepl("VarCorr|Std\\.Dev", out)))
  expect_true(any(grepl("ICC", out)))
})


# ---------------------------------------------------------------------------- #
# C: plot.brsmm — which=7 (RE Q-Q) and which=8 (caterpillar)
# ---------------------------------------------------------------------------- #

test_that("plot.brsmm which=7 (RE Q-Q) runs without error", {
  fit <- .make_brsmm_fit()
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(fit, which = 7L))
})

test_that("plot.brsmm which=8 (RE caterpillar) runs without error", {
  fit <- .make_brsmm_fit()
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(fit, which = 8L))
})


# ---------------------------------------------------------------------------- #
# D: autoplot.brsmm — title/xlab/ylab, theme, shrinkage, type='all'
# ---------------------------------------------------------------------------- #

test_that("autoplot.brsmm title/xlab/ylab override labels", {
  skip_if_not_installed("ggplot2")
  fit <- .make_brsmm_fit()
  p <- autoplot.brsmm(fit,
    type = "calibration",
    title = "MY_TITLE", xlab = "MY_X", ylab = "MY_Y"
  )
  expect_s3_class(p, "ggplot")
  lbs <- p$labels
  expect_equal(lbs$title, "MY_TITLE")
  expect_equal(lbs$x, "MY_X")
  expect_equal(lbs$y, "MY_Y")
})

test_that("autoplot.brsmm theme arg accepts theme object and function", {
  skip_if_not_installed("ggplot2")
  fit <- .make_brsmm_fit()

  # theme object
  p1 <- autoplot.brsmm(fit,
    type = "calibration",
    theme = ggplot2::theme_bw()
  )
  expect_s3_class(p1, "ggplot")

  # theme function
  p2 <- autoplot.brsmm(fit,
    type = "calibration",
    theme = ggplot2::theme_bw
  )
  expect_s3_class(p2, "ggplot")
})

test_that("autoplot.brsmm type='shrinkage' returns a ggplot", {
  skip_if_not_installed("ggplot2")
  fit <- .make_brsmm_fit()
  p <- autoplot.brsmm(fit, type = "shrinkage")
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.brsmm type='all' returns list of ggplots", {
  skip_if_not_installed("ggplot2")
  fit <- .make_brsmm_fit()
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  result <- autoplot.brsmm(fit, type = "all")
  expect_type(result, "list")
  expect_true(length(result) >= 2L)
  expect_true(all(vapply(result, inherits, logical(1L), "ggplot")))
})

test_that("autoplot.brsmm ... forwards to ggplot2::theme()", {
  skip_if_not_installed("ggplot2")
  fit <- .make_brsmm_fit()
  # Should not error even if passing valid theme args
  p <- expect_no_error(
    autoplot.brsmm(fit,
      type = "calibration",
      legend.position = "none"
    )
  )
  expect_s3_class(p, "ggplot")
})


# ---------------------------------------------------------------------------- #
# E: autoplot.brs — title/xlab/ylab, theme, type='all'
# ---------------------------------------------------------------------------- #

test_that("autoplot.brs title/xlab/ylab override labels", {
  skip_if_not_installed("ggplot2")
  fit <- .make_brs_fit()
  p <- autoplot.brs(fit,
    type = "calibration",
    title = "MY_TITLE", xlab = "MY_X", ylab = "MY_Y"
  )
  expect_s3_class(p, "ggplot")
  lbs <- p$labels
  expect_equal(lbs$title, "MY_TITLE")
  expect_equal(lbs$x, "MY_X")
  expect_equal(lbs$y, "MY_Y")
})

test_that("autoplot.brs theme arg accepted as object and function", {
  skip_if_not_installed("ggplot2")
  fit <- .make_brs_fit()
  p1 <- autoplot.brs(fit,
    type = "score_dist",
    theme = ggplot2::theme_classic()
  )
  expect_s3_class(p1, "ggplot")
  p2 <- autoplot.brs(fit,
    type = "score_dist",
    theme = ggplot2::theme_classic
  )
  expect_s3_class(p2, "ggplot")
})

test_that("autoplot.brs type='all' returns list of ggplots", {
  skip_if_not_installed("ggplot2")
  fit <- .make_brs_fit()
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  result <- autoplot.brs(fit, type = "all")
  expect_type(result, "list")
  expect_true(length(result) >= 2L)
  expect_true(all(vapply(result, inherits, logical(1L), "ggplot")))
})

test_that("autoplot.brs ... forwards to ggplot2::theme()", {
  skip_if_not_installed("ggplot2")
  fit <- .make_brs_fit()
  p <- expect_no_error(
    autoplot.brs(fit,
      type = "calibration",
      legend.position = "none"
    )
  )
  expect_s3_class(p, "ggplot")
})

test_that("autoplot.brs type='residuals_by_delta' returns ggplot", {
  skip_if_not_installed("ggplot2")
  fit <- .make_brs_fit()
  p <- autoplot.brs(fit, type = "residuals_by_delta")
  expect_s3_class(p, "ggplot")
})
