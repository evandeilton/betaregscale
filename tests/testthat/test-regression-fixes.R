testthat::test_that("brs_repar valida dominio e comprimentos", {
  testthat::expect_error(
    brs_repar(mu = c(0.3, 0.4), phi = c(0.2), repar = 2),
    NA
  )
  testthat::expect_error(
    brs_repar(mu = c(0.3, 0.4), phi = c(0.2, 0.3, 0.4), repar = 2),
    "compatible lengths"
  )
  testthat::expect_error(
    brs_repar(mu = c(0, 0.4), phi = c(0.2, 0.3), repar = 2),
    "must lie in \\(0, 1\\)"
  )
  testthat::expect_error(
    brs_repar(mu = c(0.3, 0.4), phi = c(1.2, 0.3), repar = 2),
    "repar = 2"
  )
})

testthat::test_that(".extract_response alinha prepared data com rownames nao numericos", {
  dat <- data.frame(
    y = c(10, 20, 30),
    x = c(1, 2, 3),
    left = c(0.1, 0.2, 0.3),
    right = c(0.2, 0.3, 0.4),
    yt = c(0.15, 0.25, 0.35),
    delta = c(0L, 1L, 3L),
    stringsAsFactors = FALSE
  )
  rownames(dat) <- c("a", "b", "c")
  attr(dat, "is_prepared") <- TRUE

  mf <- stats::model.frame(y ~ x, data = dat)
  ymat <- .extract_response(mf, dat, ncuts = 100L, lim = 0.5)

  testthat::expect_equal(nrow(ymat), 3L)
  testthat::expect_false(anyNA(ymat))
  testthat::expect_equal(as.integer(ymat[, "delta"]), c(0L, 1L, 3L))
})

testthat::test_that("fitted/predict de precision retornam comprimentos corretos", {
  set.seed(42)
  n <- 80
  dat <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  sim <- brs_sim(
    formula = ~ x1 + x2,
    data = dat,
    beta = c(0.2, -0.3, 0.1),
    phi = 0.2,
    ncuts = 20L
  )
  fit <- brs(y ~ x1 + x2, data = sim, ncuts = 20L)

  testthat::expect_length(stats::fitted(fit, type = "phi"), n)
  testthat::expect_length(stats::predict(fit, type = "precision"), n)

  newd <- sim[1:7, c("y", "x1", "x2")]
  testthat::expect_length(
    stats::predict(fit, newdata = newd, type = "precision"),
    7L
  )
})

testthat::test_that("residuals rqr funciona com observacoes censuradas", {
  set.seed(123)
  n <- 60
  dat <- data.frame(x1 = stats::rnorm(n))
  sim <- brs_sim(
    formula = ~x1,
    data = dat,
    beta = c(0.1, 0.2),
    phi = 0.2,
    ncuts = 15L,
    delta = 3L
  )
  fit <- brs(y ~ x1, data = sim, ncuts = 15L)

  set.seed(999)
  rq <- stats::residuals(fit, type = "rqr")
  testthat::expect_length(rq, n)
  testthat::expect_true(all(is.finite(rq)))
})
