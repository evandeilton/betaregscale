# ============================================================================ #
# Regression tests for the 2026-09 audit findings (named, known cases)
# ============================================================================ #

.audit_prep <- function(seed = 11L, n = 200L, g = 20L) {
  set.seed(seed)
  x <- rnorm(n)
  id <- factor(rep(seq_len(g), length.out = n))
  b <- rnorm(g, sd = 0.6)
  mu <- plogis(0.3 + 0.8 * x + b[as.integer(id)])
  shp <- brs_repar(mu = mu, phi = 0.2, repar = 2L)
  y <- round(stats::rbeta(n, shp$shape1, shp$shape2) * 10)
  brs_prep(data.frame(y = y, x = x, id = id), ncuts = 10L)
}

# C1 -------------------------------------------------------------------------- #

test_that("C1: brs on subsetted prepared data keeps rows aligned", {
  p <- suppressMessages(.audit_prep())
  sub <- p[-(1:50), ]
  fit_sub <- brs(y ~ x, data = sub, ncuts = 10L)
  # Fitting the same rows after resetting rownames must give the same answer
  sub2 <- sub
  rownames(sub2) <- NULL
  fit_ref <- brs(y ~ x, data = sub2, ncuts = 10L)
  expect_equal(unname(fit_sub$Y[, "left"]), unname(fit_ref$Y[, "left"]))
  expect_equal(coef(fit_sub), coef(fit_ref), tolerance = 1e-6)
  # Named case from the audit: dropping one row (jackknife) flipped the slope
  fit_jk <- brs(y ~ x, data = p[-10, ], ncuts = 10L)
  expect_true(abs(coef(fit_jk)[["x"]] - coef(fit_ref)[["x"]]) < 0.3)
})

test_that("C1: brs_cv works on data from brs_prep", {
  p <- suppressMessages(.audit_prep())
  set.seed(1L)
  cv <- brs_cv(y ~ x, data = p, k = 3L, ncuts = 10L)
  expect_s3_class(cv, "brs_cv")
  expect_true(all(is.finite(cv$log_score)))
})

test_that("C1: brsmm is invariant to row permutation and subsetting", {
  p <- suppressMessages(.audit_prep())
  ctrl <- list(maxit = 300L)
  fit0 <- brsmm(y ~ x, random = ~ 1 | id, data = p, ncuts = 10L,
                control = ctrl)
  set.seed(3L)
  fit_perm <- brsmm(y ~ x, random = ~ 1 | id, data = p[sample(nrow(p)), ],
                    ncuts = 10L, control = ctrl)
  expect_equal(fit_perm$value, fit0$value, tolerance = 1e-4)
  expect_equal(fit_perm$random$sigma_b, fit0$random$sigma_b, tolerance = 1e-2)
  # Every-other-row subset used to error with "Grouping variable contains
  # missing values"; it must simply fit.
  fit_odd <- brsmm(y ~ x, random = ~ 1 | id, data = p[seq(1, nrow(p), 2), ],
                   ncuts = 10L, control = ctrl)
  expect_s3_class(fit_odd, "brsmm")
})

# C3 -------------------------------------------------------------------------- #

test_that("C3: compiled brsmm likelihood rejects inconsistent dimensions", {
  p <- suppressMessages(.audit_prep(n = 40L, g = 4L))
  X <- cbind(1, p$x); Z <- matrix(1, nrow(p), 1); Xr <- matrix(1, nrow(p), 1)
  par <- c(0, 0, 0, log(0.5))
  args <- list(par, X, Z, Xr, p$left, p$right, p$yt, as.integer(p$delta),
               as.integer(p$id), 0L, 0L, 2L, 0L, 11L)
  expect_true(is.finite(do.call(betaregscale:::.brsmm_loglik_eigen, args)))
  bad <- args; bad[[9]] <- as.integer(p$id)[-1]
  expect_error(do.call(betaregscale:::.brsmm_loglik_eigen, bad), "rows")
  bad <- args; bad[[9]][1L] <- 0L
  expect_error(do.call(betaregscale:::.brsmm_loglik_eigen, bad), ">= 1")
  bad <- args; bad[[4]] <- Xr[-1, , drop = FALSE]
  expect_error(do.call(betaregscale:::.brsmm_loglik_eigen, bad), "rows")
})

test_that("C3: NA in a random-slope variable gives a clear error", {
  p <- suppressMessages(.audit_prep())
  p$z <- rnorm(nrow(p)); p$z[nrow(p)] <- NA
  expect_error(
    brsmm(y ~ x, random = ~ z | id, data = p, ncuts = 10L),
    "missing values"
  )
})

# C4 -------------------------------------------------------------------------- #

test_that("C4: Mode 3 (left/right only) observations are kept in the fit", {
  set.seed(5L)
  n <- 120L
  x <- rnorm(n)
  y <- round(plogis(0.2 + 0.5 * x) * 10 + rnorm(n))
  y <- pmin(pmax(y, 0), 10)
  d <- data.frame(y = y, x = x, left = NA_real_, right = NA_real_)
  # 40 rows become Mode 3: only a right bound is known (left-censored)
  m3 <- 1:40
  d$y[m3] <- NA
  d$right[m3] <- pmax(y[m3], 1)
  p <- suppressMessages(brs_prep(d, ncuts = 10L))
  expect_equal(sum(p$delta == 1L), 40L)
  expect_false(anyNA(p$y))
  fit <- brs(y ~ x, data = p, ncuts = 10L)
  expect_equal(nobs(fit), n)
  expect_equal(sum(fit$Y[, "delta"] == 1L), 40L)
})
