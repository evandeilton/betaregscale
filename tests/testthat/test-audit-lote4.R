# ============================================================================ #
# Lote 4 (2026-09 audit): interval = c("mid", "right", "left")
# Named cases L4-*. Identities use 1e-10 (or 1e-8 across fits), recovery 3 SE.
# ============================================================================ #

.l4_eps <- 1e-5

# Simulated scores (repar 2, K = 10 by default) with the given direction.
.l4_sim <- function(interval, n = 300L, seed = 41L, K = 10L, beta = c(0.3, 0.8),
                    phi = -1.5) {
  set.seed(seed)
  d <- data.frame(x = runif(n, -1, 1))
  brs_sim(y ~ x, data = d, beta = beta, phi = phi, ncuts = K, interval = interval)
}

# Raw (non-prepared) scores of a simulated data set.
.l4_raw <- function(s) data.frame(y = s$y, x = s$x)

# L4-1: cell mapping ----------------------------------------------------------

test_that("L4-1a: brs_check(0:10, ncuts = 10) cells per interval (snapshot for mid)", {
  s <- 0:10
  r <- brs_check(s, ncuts = 10, interval = "right")
  expect_equal(unname(r[, "left"]), pmax(s / 11, .l4_eps), tolerance = 1e-10)
  expect_equal(unname(r[, "right"]), pmin((s + 1) / 11, 1 - .l4_eps), tolerance = 1e-10)
  expect_equal(unname(r[, "yt"]), (s + 0.5) / 11, tolerance = 1e-10)
  expect_identical(unname(r[, "delta"]), c(1, rep(3, 9), 2))
  l <- brs_check(s, ncuts = 10, interval = "left")
  expect_identical(l, r)
  # "mid" is today's behaviour, unchanged
  m <- brs_check(s, ncuts = 10)
  expect_identical(m, brs_check(s, ncuts = 10, interval = "mid"))
  expect_equal(unname(m[, "left"]), pmax((s - 0.5) / 10, .l4_eps), tolerance = 1e-10)
  expect_equal(unname(m[, "right"]), pmin((s + 0.5) / 10, 1 - .l4_eps), tolerance = 1e-10)
  expect_equal(unname(m[, "yt"]), pmin(pmax(s / 10, .l4_eps), 1 - .l4_eps), tolerance = 1e-10)
  expect_identical(unname(m[, "delta"]), c(1, rep(3, 9), 2))
  # snapshot of a mid row: score 3 of 10 -> [0.25, 0.35], yt 0.3
  expect_equal(unname(m[4, c("left", "right", "yt")]), c(0.25, 0.35, 0.3), tolerance = 1e-12)
})

test_that("L4-1b: edge grids: K = 1, K = 2, one score, one border, mixed proportions", {
  # K = 1: two border cells [0, 1/2], [1/2, 1] in every mode
  for (iv in c("mid", "right", "left")) {
    r <- brs_check(c(0, 1), ncuts = 1, interval = iv)
    expect_identical(unname(r[, "delta"]), c(1, 2))
    expect_equal(unname(r[, "right"])[1], 0.5, tolerance = 1e-10)
    expect_equal(unname(r[, "left"])[2], 0.5, tolerance = 1e-10)
  }
  # K = 2
  r2 <- brs_check(0:2, ncuts = 2, interval = "right")
  expect_equal(unname(r2[, "left"]), c(.l4_eps, 1 / 3, 2 / 3), tolerance = 1e-10)
  expect_equal(unname(r2[, "right"]), c(1 / 3, 2 / 3, 1 - .l4_eps), tolerance = 1e-10)
  m2 <- brs_check(0:2, ncuts = 2)
  expect_equal(unname(m2[, "left"]), c(.l4_eps, 0.25, 0.75), tolerance = 1e-10)
  expect_equal(unname(m2[, "right"]), c(0.25, 0.75, 1 - .l4_eps), tolerance = 1e-10)
  # single score value and all scores at one border
  one <- brs_check(rep(5, 4), ncuts = 10, interval = "right")
  expect_true(all(one[, "delta"] == 3) && all(one[, "left"] == 5 / 11))
  b0 <- brs_check(rep(0, 4), ncuts = 10, interval = "left")
  expect_true(all(b0[, "delta"] == 1) && all(b0[, "right"] == 1 / 11))
  bK <- brs_check(rep(10, 4), ncuts = 10, interval = "right")
  expect_true(all(bK[, "delta"] == 2) && all(bK[, "left"] == 10 / 11))
  # a response entirely in (0, 1) is exact in every mode
  for (iv in c("mid", "right", "left")) {
    u <- suppressMessages(brs_check(c(0.2, 0.7), ncuts = 10, interval = iv))
    expect_true(all(u[, "delta"] == 0) && all(u[, "left"] == c(0.2, 0.7)))
  }
  # proportions mixed with scores: today's rule (fractional score) is kept
  mx <- brs_check(c(0.3, 5), ncuts = 10)
  expect_identical(unname(mx[, "delta"]), c(3, 3))
  expect_equal(unname(mx[1, "right"]), 0.08, tolerance = 1e-10)
  mxr <- brs_check(c(0.3, 5), ncuts = 10, interval = "right")
  expect_equal(unname(mxr[1, c("left", "right")]), c(0.3, 1.3) / 11, tolerance = 1e-10)
  # forced delta = 0 makes a proportion exact
  mx0 <- brs_check(c(0.3, 5), ncuts = 10, delta = c(0L, 3L), interval = "right")
  expect_equal(unname(mx0[1, c("left", "right", "yt")]), rep(0.3, 3), tolerance = 1e-10)
})

test_that("L4-1c: forced delta under right/left keeps the cell of the score", {
  r <- brs_check(rep(3, 4), ncuts = 10, delta = c(0L, 1L, 2L, 3L), interval = "right")
  expect_equal(unname(r[1, c("left", "right", "yt")]), rep(3.5 / 11, 3), tolerance = 1e-10)
  expect_equal(unname(r[2, c("left", "right")]), c(.l4_eps, 4 / 11), tolerance = 1e-10)
  expect_equal(unname(r[3, c("left", "right")]), c(3 / 11, 1 - .l4_eps), tolerance = 1e-10)
  expect_equal(unname(r[4, c("left", "right")]), c(3 / 11, 4 / 11), tolerance = 1e-10)
  expect_identical(unname(r[, "delta"]), c(0, 1, 2, 3))
  # mid, forced delta: unchanged formulas u = (y + 0.5)/K, l = (y - 0.5)/K
  m <- brs_check(c(30, 60), ncuts = 100, delta = c(1L, 2L))
  expect_equal(unname(m[1, c("left", "right")]), c(.l4_eps, 0.305), tolerance = 1e-10)
  expect_equal(unname(m[2, c("left", "right")]), c(0.595, 1 - .l4_eps), tolerance = 1e-10)
})

test_that("L4-1d: D2 - delta = 3 with left == right is rejected", {
  # a border cell squeezed to [eps, eps] by the clamp is reported as such (item 10)
  expect_error(
    suppressWarnings(brs_check(c(0, 5), ncuts = 10, lim = 1e-9, delta = c(3L, 3L))),
    "narrower than the 1e-5 border clamp"
  )
  expect_error(
    suppressMessages(brs_prep(data.frame(y = 50, left = 50, right = 50), ncuts = 100)),
    "delta = 3 with left >= right"
  )
  # the same point declared exact is fine
  ok <- suppressMessages(brs_prep(
    data.frame(y = 50, left = 50, right = 50, delta = 0), ncuts = 100))
  expect_identical(ok$delta, 0L)
  expect_equal(ok$yt, 0.5, tolerance = 1e-10)
})

test_that("L4-1f: scores outside 0..K stop in brs_check() and in the fit (not via D2)", {
  for (iv in c("mid", "right", "left")) {
    expect_error(brs_check(c(5, 11), ncuts = 10, interval = iv),
                 "Maximum response \\(11\\) exceeds ncuts \\(10\\)")
    expect_error(brs_check(c(5, 11), ncuts = 10, delta = c(3L, 3L), interval = iv),
                 "exceeds ncuts")
    expect_error(brs_check(c(-1, 5), ncuts = 10, interval = iv), "non-negative")
  }
  d <- data.frame(y = c(0, 4, 12, 7), x = 1:4)
  expect_error(brs(y ~ x, data = d, ncuts = 10), "Increase `ncuts` to at least 12")
  expect_error(brs(y ~ x, data = d, ncuts = 10, interval = "right"), "exceeds ncuts")
  # brs_prep() applies the same rule
  expect_error(brs_prep(d, ncuts = 10), "must be >= the maximum observed value")
})

test_that("L4-1e: brs_prep() Modes 1-4 under right/left and the interval attribute", {
  p <- suppressMessages(brs_prep(data.frame(y = c(0, 3, 10), x = 1:3), ncuts = 10,
                                 interval = "right"))
  expect_identical(attr(p, "interval"), "right")
  expect_equal(p$left, pmax(c(0, 3, 10) / 11, .l4_eps), tolerance = 1e-10)
  expect_equal(p$right, pmin(c(1, 4, 11) / 11, 1 - .l4_eps), tolerance = 1e-10)
  expect_identical(p$delta, c(1L, 3L, 2L))
  # Mode 2 (explicit delta) matches brs_check with forced delta
  p2 <- suppressMessages(suppressWarnings(brs_prep(
    data.frame(y = c(3, 3, 3), delta = c(1, 2, 3)), ncuts = 10, interval = "left")))
  chk <- brs_check(c(3, 3, 3), ncuts = 10, delta = c(1L, 2L, 3L), interval = "left")
  expect_equal(p2$left, unname(chk[, "left"]), tolerance = 1e-10)
  expect_equal(p2$right, unname(chk[, "right"]), tolerance = 1e-10)
  # Modes 3/4: analyst endpoints are latent scores: L/K, L/(K+1), (L+1)/(K+1)
  p4 <- suppressMessages(brs_prep(data.frame(left = 30, right = 45), ncuts = 100,
                                  interval = "right"))
  expect_equal(c(p4$left, p4$right), c(30, 45) / 101, tolerance = 1e-10)
  p4l <- suppressMessages(brs_prep(data.frame(left = 30, right = 45), ncuts = 100,
                                   interval = "left"))
  expect_equal(c(p4l$left, p4l$right), c(31, 46) / 101, tolerance = 1e-10)
  p4m <- suppressMessages(brs_prep(data.frame(left = 30, right = 45), ncuts = 100))
  expect_equal(c(p4m$left, p4m$right), c(0.30, 0.45), tolerance = 1e-10)
  expect_identical(attr(p4m, "interval"), "mid")
  # the dissertation's interval of score 50 (m, r, l) gives the cell of score 50
  bounds <- list(mid = c(49.5, 50.5), right = c(50, 51), left = c(49, 50))
  for (iv in names(bounds)) {
    pa <- suppressMessages(brs_prep(
      data.frame(y = 50, left = bounds[[iv]][1], right = bounds[[iv]][2]),
      ncuts = 100, interval = iv))
    ps <- suppressMessages(brs_prep(data.frame(y = 50), ncuts = 100, interval = iv))
    expect_equal(c(pa$left, pa$right), c(ps$left, ps$right), tolerance = 1e-12,
                 label = iv)
  }
  # Mode 3 fill of `y`: latent score of the cell centre (read back as score 50)
  p3 <- suppressMessages(brs_prep(data.frame(left = c(NA, 49), right = c(5, 50)),
                                  ncuts = 100, interval = "left"))
  expect_equal(p3$y[2], 49.5, tolerance = 1e-10)
  expect_equal(betaregscale:::.brs_observed_scores(p3$y, 100L, "left")[2], 50)
})

# L4-2: score probabilities sum to 1 ------------------------------------------

test_that("L4-2: rowSums(brs_predict_scoreprob()) == 1 in the three modes (lim = 0.5)", {
  for (K in c(10L, 2L)) {
    for (iv in c("mid", "right", "left")) {
      s <- .l4_sim(iv, n = 120L, K = K)
      f <- brs(y ~ x, data = s)
      expect_identical(f$interval, iv)
      P <- brs_predict_scoreprob(f)
      expect_equal(unname(rowSums(P)), rep(1, nrow(P)), tolerance = 1e-8,
                   label = paste("K", K, iv))
      Pn <- brs_predict_scoreprob(f, newdata = data.frame(x = c(-1, 0, 1)))
      expect_equal(unname(rowSums(Pn)), rep(1, 3), tolerance = 1e-8)
      expect_equal(ncol(P), K + 1L)
    }
  }
})

# L4-3: parameter recovery per mode -------------------------------------------

test_that("L4-3: brs_sim + brs recover beta and phi within 3 SE in each mode (n = 2000)", {
  truth <- c(0.3, 0.8, -1.5)
  for (iv in c("mid", "right", "left")) {
    s <- .l4_sim(iv, n = 2000L, seed = 2026L)
    f <- brs(y ~ x, data = s)
    se <- sqrt(diag(vcov(f)))
    z <- (coef(f) - truth) / se
    expect_true(all(abs(z) < 3), info = paste(iv, ":", paste(round(z, 2), collapse = " ")))
    expect_identical(f$convergence, 0L)
  }
})

# L4-4: right and left are the same likelihood -------------------------------

test_that("L4-4: right and left give identical fits; latent scores differ by one", {
  s <- .l4_sim("right", n = 300L)
  raw <- .l4_raw(s)
  fr <- brs(y ~ x, data = raw, ncuts = 10, interval = "right")
  fl <- brs(y ~ x, data = raw, ncuts = 10, interval = "left")
  expect_equal(as.numeric(logLik(fr)), as.numeric(logLik(fl)), tolerance = 1e-8)
  expect_equal(coef(fr), coef(fl), tolerance = 1e-8)
  expect_identical(fr$Y, fl$Y)
  sr <- predict(fr, type = "score")
  sl <- predict(fl, type = "score")
  expect_equal(unname(sl), unname(sr) - 1, tolerance = 1e-12)
  expect_equal(unname(sr), unname(11 * fitted(fr)), tolerance = 1e-12)
  expect_equal(predict(fr, type = "expected_score"), predict(fl, type = "expected_score"),
               tolerance = 1e-12)
  # mid: latent score K * E[Y]
  fm <- brs(y ~ x, data = raw, ncuts = 10)
  expect_equal(unname(predict(fm, type = "score")), unname(10 * fitted(fm)), tolerance = 1e-12)
  es <- predict(fm, type = "expected_score", newdata = data.frame(x = c(-1, 1)))
  expect_true(all(es >= 0 & es <= 10) && es[2] > es[1])
  # expected score = sum s P(s) from the score probabilities
  P <- brs_predict_scoreprob(fr)
  expect_equal(unname(predict(fr, type = "expected_score")), unname(as.numeric(P %*% 0:10)),
               tolerance = 1e-12)
})

# L4-5: lim validation ---------------------------------------------------------

test_that("L4-5: lim > 0.5 errors, lim < 0.5 warns, non-default lim is ignored under right/left", {
  s <- .l4_sim("mid", n = 60L)
  raw <- .l4_raw(s)
  expect_error(brs(y ~ x, data = raw, ncuts = 10, lim = 1), "\\(0, 0.5\\]")
  expect_error(brs(y ~ x, data = raw, ncuts = 10, lim = 0), "\\(0, 0.5\\]")
  expect_error(brs_check(0:10, ncuts = 10, lim = 0.6), "\\(0, 0.5\\]")
  expect_warning(f25 <- brs(y ~ x, data = raw, ncuts = 10, lim = 0.25), "cover only 0.5")
  expect_equal(f25$lim, 0.25)
  expect_lt(max(rowSums(brs_predict_scoreprob(f25))), 0.9)
  expect_warning(brs_check(0:10, ncuts = 10, lim = 0.25), "probabilities do not sum to 1")
  expect_warning(brs_check(0:10, ncuts = 10, lim = 0.3, interval = "right"), "ignored")
  expect_warning(fr <- brs(y ~ x, data = raw, ncuts = 10, lim = 0.3, interval = "right"),
                 "ignored")
  # the ignored lim does not change the right cells
  expect_identical(fr$Y, brs(y ~ x, data = raw, ncuts = 10, interval = "right")$Y)
  expect_warning(
    brs_sim(y ~ x, data = data.frame(x = rnorm(20)), beta = c(0, 0.5), phi = -1,
            ncuts = 10, lim = 0.25),
    "rounds to the nearest score"
  )
  # a stored lim < 0.5 warns once, in brs_prep(), not again in brs()
  expect_warning(p25 <- suppressMessages(
    brs_prep(data.frame(y = raw$y, x = raw$x), ncuts = 10, lim = 0.25)), "cover only")
  expect_silent(brs(y ~ x, data = p25))
})

# L4-6: anova refuses different coarsenings -----------------------------------

test_that("L4-6: anova() refuses fits with different interval, ncuts or lim", {
  s <- .l4_sim("mid", n = 150L)
  raw <- .l4_raw(s)
  f_mid <- brs(y ~ x, data = raw, ncuts = 10)
  f_right <- brs(y ~ x, data = raw, ncuts = 10, interval = "right")
  expect_error(anova(f_mid, f_right), "same `interval`")
  f_k20 <- suppressWarnings(brs(y ~ x, data = raw, ncuts = 20))
  expect_error(anova(f_mid, f_k20), "same `ncuts`")
  f_l25 <- suppressWarnings(brs(y ~ x, data = raw, ncuts = 10, lim = 0.25))
  expect_error(anova(f_mid, f_l25), "same `lim`")
  # nested fits with the same coarsening still compare
  f0 <- brs(y ~ 1, data = raw, ncuts = 10, interval = "right")
  a <- anova(f0, f_right)
  expect_s3_class(a, "anova")
  expect_true(a$Chisq[2] >= 0)
  # lim is ignored under right/left, so it does not block the comparison there
  f0_l3 <- suppressWarnings(brs(y ~ 1, data = raw, ncuts = 10, lim = 0.3,
                                interval = "right"))
  expect_equal(as.numeric(logLik(f0_l3)), as.numeric(logLik(f0)), tolerance = 1e-10)
  expect_s3_class(anova(f0_l3, f_right), "anova")
})

# L4-7: attributes from brs_prep()/brs_sim() ----------------------------------

test_that("L4-7: prepared interval wins over an explicit different value, with a warning", {
  s <- .l4_sim("mid", n = 100L)
  p <- suppressMessages(brs_prep(.l4_raw(s), ncuts = 10, interval = "right"))
  expect_warning(f <- brs(y ~ x, data = p, interval = "mid"), "`interval = mid` differs")
  expect_identical(f$interval, "right")
  expect_identical(brs(y ~ x, data = p)$interval, "right")
  expect_silent(brs(y ~ x, data = p, interval = "right"))
  # brs_sim() attaches interval; brs()/brsmm()/brs_cv()/bootstrap reuse it
  sl <- .l4_sim("left", n = 120L)
  expect_identical(attr(sl, "interval"), "left")
  fl <- brs(y ~ x, data = sl)
  expect_identical(fl$interval, "left")
  expect_identical(unname(fl$Y[, "left"]),
                   unname(brs_check(sl$y, ncuts = 10, interval = "left")[, "left"]))
  sl$id <- factor(rep(1:12, each = 10))
  fm <- brsmm(y ~ x, random = ~ 1 | id, data = sl, control = list(maxit = 150L))
  expect_identical(fm$interval, "left")
  expect_equal(unname(rowSums(betaregscale:::.brs_score_prob_matrix(
    fm$fitted_mu, fm$fitted_phi, 2L, 10L, 0.5, 0:10, "left"))), rep(1, nrow(sl)),
    tolerance = 1e-8)
  expect_error(brs(y ~ x, data = .l4_raw(s), interval = "middle"), "must be one of")
  # partial matching, as match.arg() in brs_check()/brs_prep()
  expect_identical(brs(y ~ x, data = .l4_raw(s), ncuts = 10, interval = "r")$interval,
                   "right")
})

# L4-8: bootstrap and cross-validation under "right" --------------------------

test_that("L4-8: brs_bootstrap and brs_cv work under interval = 'right'", {
  s <- .l4_sim("right", n = 150L, seed = 7L)
  f <- brs(y ~ x, data = s)
  set.seed(1L)
  bt <- brs_bootstrap(f, R = 20L, keep_draws = TRUE)
  expect_identical(attr(bt, "n_success"), 20L)
  draws <- attr(bt, "boot_draws")
  z <- (colMeans(draws) - f$par) / bt$se_boot
  expect_true(all(abs(z) < 3), info = paste(round(z, 2), collapse = " "))
  # cv metrics with train == test: mean per-observation in-sample loglik
  m <- betaregscale:::.brs_cv_metrics(f, f$data)
  expect_equal(m$log_score, as.numeric(logLik(f)) / nobs(f), tolerance = 1e-8)
  set.seed(2L)
  cv <- brs_cv(y ~ x, data = s, k = 3L)
  expect_true(all(cv$converged) && all(is.finite(cv$log_score)))
  # raw scores with interval forwarded through ...
  cv2 <- brs_cv(y ~ x, data = .l4_raw(s), k = 3L, ncuts = 10, interval = "right")
  expect_true(all(is.finite(cv2$log_score)))
})

# L4-9: DGP of brs_sim matches the score probabilities -----------------------

test_that("L4-9: empirical score frequencies of brs_sim match mean_i P_i(s) (chi-square)", {
  for (iv in c("mid", "right")) {
    set.seed(99L)
    n <- 2000L
    K <- 10L
    x <- runif(n, -1, 1)
    s <- brs_sim(y ~ x, data = data.frame(x), beta = c(0.3, 0.8), phi = -1.5,
                 ncuts = K, interval = iv)
    P <- betaregscale:::.brs_score_prob_matrix(
      plogis(0.3 + 0.8 * x), plogis(-1.5), 2L, K, 0.5, 0:K, iv)
    expected <- colSums(P)
    observed <- tabulate(s$y + 1L, K + 1L)
    chisq <- sum((observed - expected)^2 / expected)
    pval <- pchisq(chisq, df = K, lower.tail = FALSE)
    expect_gt(pval, 1e-3, label = paste(iv, "chisq", round(chisq, 2)))
    expect_true(all(expected > 5))
    expect_equal(sum(expected), n, tolerance = 1e-8)
  }
})

# L4-10: brsmm under "right": AGHQ vs brute-force integration -----------------

.l4_marginal_ll <- function(fit) {
  X <- fit$model_matrices$X
  Z <- fit$model_matrices$Z
  Xr <- fit$model_matrices$Xr
  beta <- fit$par[seq_len(fit$p)]
  gamma <- fit$par[fit$p + seq_len(fit$q)]
  sd_b <- as.numeric(fit$random$sd_b[1L])
  eta_fixed <- as.numeric(X %*% beta)
  phi <- betaregscale:::.clamp_phi_by_repar(
    betaregscale:::apply_inv_link(as.numeric(Z %*% gamma), fit$link_phi), fit$repar)
  Y <- fit$Y
  g <- fit$group_index
  total <- 0
  for (k in seq_len(fit$ngroups)) {
    idx <- which(g == k)
    f <- function(b) {
      vapply(b, function(bb) {
        mu <- betaregscale:::.clamp_mu_by_repar(
          betaregscale:::apply_inv_link(eta_fixed[idx] + Xr[idx, 1L] * bb, fit$link),
          fit$repar)
        sh <- brs_repar(mu, phi[idx], repar = fit$repar)
        exp(sum(betaregscale:::.brs_obs_loglik(
          Y[idx, "delta"], Y[idx, "left"], Y[idx, "right"], Y[idx, "yt"],
          sh$shape1, sh$shape2)) + dnorm(bb, 0, sd_b, log = TRUE))
      }, numeric(1))
    }
    # finite range: integrate() over (-Inf, Inf) was only accurate to ~1e-7 relative
    total <- total + log(integrate(f, -12 * sd_b, 12 * sd_b, rel.tol = 1e-12,
                                   subdivisions = 1000L)$value)
  }
  total
}

test_that("L4-10: brsmm AGHQ log-likelihood under 'right' equals brute-force integration", {
  skip_on_cran()
  set.seed(5L)
  g <- 12L
  ni <- 10L
  n <- g * ni
  id <- factor(rep(seq_len(g), each = ni))
  x <- runif(n, -1, 1)
  b <- rnorm(g, sd = 0.5)
  mu <- plogis(0.3 + 0.8 * x + b[as.integer(id)])
  shp <- brs_repar(mu, 0.2, repar = 2L)
  y <- betaregscale:::.brs_score_from_unit(rbeta(n, shp$shape1, shp$shape2), 10L, "right")
  d <- data.frame(y = y, x = x, id = id)
  f <- brsmm(y ~ x, random = ~ 1 | id, data = d, ncuts = 10, interval = "right",
             int_method = "aghq", n_points = 25L, control = list(maxit = 300L))
  expect_identical(f$interval, "right")
  ll_bf <- .l4_marginal_ll(f)
  expect_lt(abs(f$value - ll_bf), 1e-6)
  expect_true(all(f$Y[, "delta"] %in% 1:3))
})

# L4-11: observed scores and helpers ------------------------------------------

test_that("L4-11: .brs_observed_scores and .brs_score_from_unit follow the interval", {
  y <- c(0.3, 0.36, 0.04, 0.96)
  expect_equal(betaregscale:::.brs_observed_scores(y, 10L, "mid"), c(3, 4, 0, 10))
  expect_equal(betaregscale:::.brs_observed_scores(y, 10L, "right"), c(3, 3, 0, 10))
  expect_equal(betaregscale:::.brs_observed_scores(y, 10L), c(3, 4, 0, 10))
  expect_equal(betaregscale:::.brs_score_from_unit(c(0, 0.999999, 1), 10L, "right"), c(0, 10, 10))
  expect_equal(betaregscale:::.brs_score_from_unit(c(0, 0.96, 1), 10L, "mid"), c(0, 10, 10))
  expect_equal(betaregscale:::.brs_latent_score(0.5, 10L, "left"), 4.5)
  expect_equal(betaregscale:::.brs_latent_score(0.5, 10L, "right"), 5.5)
  expect_equal(betaregscale:::.brs_latent_score(0.5, 10L, "mid"), 5)
  expect_identical(betaregscale:::.brs_interval_of(list()), "mid")
  # integer scores in {0, 1} are scores, not unit values (single score 1, K = 10)
  expect_equal(betaregscale:::.brs_observed_scores(c(1, 1, 1), 10L), c(1, 1, 1))
  expect_equal(betaregscale:::.brs_observed_scores(c(0, 1, 0), 10L, "right"), c(0, 1, 0))
  # latent values: the cell containing them (round / floor / ceiling)
  expect_equal(betaregscale:::.brs_observed_scores(c(37.5, 2), 100L, "right"), c(37, 2))
  expect_equal(betaregscale:::.brs_observed_scores(c(37.5, 2), 100L, "left"), c(38, 2))
  # .brs_unit_from_latent inverts .brs_latent_score
  for (iv in c("mid", "right", "left")) {
    expect_equal(betaregscale:::.brs_unit_from_latent(
      betaregscale:::.brs_latent_score(0.37, 10L, iv), 10L, iv), 0.37, tolerance = 1e-12)
  }
})

# L4-12: brs_sim with forced delta under right/left --------------------------

test_that("L4-12: brs_sim forced delta under 'right' keeps the score cells", {
  set.seed(3L)
  d <- data.frame(x = runif(50, -1, 1))
  s3 <- brs_sim(y ~ x, data = d, beta = c(0.3, 0.8), phi = -1.5, ncuts = 10,
                interval = "right", delta = 3)
  expect_true(all(s3$delta == 3L) && all(s3$y >= 1 & s3$y <= 9))
  expect_equal(s3$left, s3$y / 11, tolerance = 1e-10)
  expect_equal(s3$right, (s3$y + 1) / 11, tolerance = 1e-10)
  s1 <- brs_sim(y ~ x, data = d, beta = c(0.3, 0.8), phi = -1.5, ncuts = 10,
                interval = "left", delta = 1)
  expect_true(all(s1$delta == 1L) && all(s1$left == .l4_eps))
  expect_equal(s1$right, pmin((s1$y + 1) / 11, 1 - .l4_eps), tolerance = 1e-10)
  s0 <- brs_sim(y ~ x, data = d, beta = c(0.3, 0.8), phi = -1.5, ncuts = 10,
                interval = "right", delta = 0)
  expect_true(all(s0$delta == 0L) && all(s0$left == s0$right))
  expect_identical(attr(s0, "interval"), "right")
  # every simulated data set fits under its own attributes
  expect_s3_class(brs(y ~ x, data = s3), "brs")
})

# L4-13: gradient on right/left data ------------------------------------------

test_that("L4-13: C++ gradient matches numDeriv on right-direction data", {
  s <- .l4_sim("right", n = 150L, seed = 8L)
  X <- cbind(1, s$x)
  par <- c(0.3, 0.8, -1.5) + c(0.05, -0.05, 0.1)
  fn <- function(p) betaregscale:::.brs_loglik_fixed_cpp(
    p, X, s$left, s$right, s$yt, s$delta, 0L, 0L, 2L)
  gc <- betaregscale:::.brs_grad_fixed_cpp(par, X, s$left, s$right, s$yt, s$delta, 0L, 0L, 2L)
  gn <- numDeriv::grad(fn, par)
  expect_lt(max(abs(gc - gn)) / max(abs(gn)), 1e-6)
})

# ============================================================================ #
# Validation round (L4-14 .. L4-23): consolidated validator fixes
# ============================================================================ #

test_that("L4-14: brs_prep checks analyst bounds on the latent range of the direction", {
  pm <- function(d, iv) suppressMessages(brs_prep(d, ncuts = 10, interval = iv))
  # the top cell of score 10 under right ([10, 11]) is accepted
  p <- pm(data.frame(left = 10, right = 11), "right")
  expect_identical(p$delta, 2L)
  # [-1, 0] under right: range error first, not the D2 "Use delta = 0" message
  expect_error(pm(data.frame(left = -1, right = 0), "right"),
               "outside the latent scale \\[0, 11\\] for interval = 'right'")
  expect_error(pm(data.frame(left = 11, right = NA_real_), "mid"),
               "outside the latent scale \\[-0.5, 10.5\\] for interval = 'mid'")
  expect_error(pm(data.frame(left = NA_real_, right = 11), "left"),
               "Column 'right'.*\\[-1, 10\\] for interval = 'left'")
  # bounds at the ends of each range are accepted
  expect_silent(pm(data.frame(left = -0.5, right = 10.5), "mid"))
  expect_silent(pm(data.frame(left = -1, right = 10), "left"))
  # scores keep the 0..K rule
  expect_error(pm(data.frame(y = 11), "right"), "must be >= the maximum observed value")
})

test_that("L4-15: prepared data without an interval attribute count as mid cells", {
  p <- suppressMessages(brs_prep(.l4_raw(.l4_sim("mid", n = 60L)), ncuts = 10))
  attr(p, "interval") <- NULL
  expect_warning(f <- brs(y ~ x, data = p, interval = "right"), "`interval = right` differs")
  expect_identical(f$interval, "mid")
  f0 <- brs(y ~ x, data = p)
  expect_identical(f0$interval, "mid")
  expect_equal(predict(f, type = "expected_score"), predict(f0, type = "expected_score"),
               tolerance = 1e-12)
})

test_that("L4-16: analyst intervals reaching 0 or 1 become one-sided censoring", {
  bounds <- list(mid = c(-0.5, 0.5, 9.5, 10.5), right = c(0, 1, 10, 11),
                 left = c(-1, 0, 9, 10))
  for (iv in names(bounds)) {
    b <- bounds[[iv]]
    pa <- suppressMessages(brs_prep(data.frame(left = b[c(1, 3)], right = b[c(2, 4)]),
                                    ncuts = 10, interval = iv))
    ps <- suppressMessages(brs_prep(data.frame(y = c(0, 10)), ncuts = 10, interval = iv))
    expect_identical(pa$delta, c(1L, 2L), label = iv)
    expect_equal(c(pa$left, pa$right), c(ps$left, ps$right), tolerance = 1e-12)
  }
  # an analyst delta is kept; an interval covering (0, 1) stays delta = 3
  keep <- suppressMessages(brs_prep(
    data.frame(left = -0.5, right = 0.5, delta = 3), ncuts = 10))
  expect_identical(keep$delta, 3L)
  whole <- suppressMessages(brs_prep(data.frame(left = -0.5, right = 10.5), ncuts = 10))
  expect_identical(whole$delta, 3L)
})

test_that("L4-17: a proportion with a forced delta keeps the score-based yt of ba5abf3", {
  p <- suppressWarnings(suppressMessages(brs_prep(
    data.frame(y = rep(0.3, 4), delta = 0:3), ncuts = 100)))
  expect_equal(p$yt, c(0.3, 0.003, 0.003, 0.003), tolerance = 1e-12)
  expect_equal(p$right[c(2, 4)], c(0.008, 0.008), tolerance = 1e-12)
})

test_that("L4-18: vectorised brs_prep agrees with brs_check for every direction and delta", {
  set.seed(18L)
  y <- sample(0:10, 200, TRUE)
  dl <- sample(0:3, 200, TRUE)
  dl[dl == 3L & y %in% c(0, 10)] <- 0L          # keep the D2 rule out of the way
  for (iv in c("mid", "right", "left")) {
    p <- suppressWarnings(suppressMessages(brs_prep(data.frame(y = y, delta = dl),
                                                    ncuts = 10, interval = iv)))
    chk <- suppressWarnings(brs_check(y, ncuts = 10, delta = dl, interval = iv))
    expect_equal(p$left, unname(chk[, "left"]), tolerance = 1e-12, label = iv)
    expect_equal(p$right, unname(chk[, "right"]), tolerance = 1e-12, label = iv)
    expect_identical(p$delta, as.integer(chk[, "delta"]))
  }
})

test_that("L4-19: lim within rounding of 0.5 counts as 0.5; messages print full digits", {
  expect_silent(brs_check(0:10, ncuts = 10, lim = 0.7 - 0.2))
  expect_silent(brs_check(0:10, ncuts = 10, lim = 0.5 + 1e-12, interval = "right"))
  expect_error(brs_check(0:10, ncuts = 10, lim = 0.5 + 1e-6), "\\(0, 0.5\\]")
  expect_warning(brs_check(0:10, ncuts = 10, lim = 0.123456789), "lim = 0.123456789")
})

test_that("L4-20: internal entry points match interval by name", {
  s <- .l4_sim("mid", n = 60L)
  raw <- .l4_raw(s)
  par <- c(0.3, 0.8, -1.5)
  ll_mid <- betaregscale:::brs_loglik(par, y ~ x, raw, ncuts = 10, interval = "mid")
  expect_equal(betaregscale:::brs_loglik(par, y ~ x, raw, ncuts = 10, interval = "m"),
               ll_mid, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(
    betaregscale:::brs_loglik(par, y ~ x, raw, ncuts = 10, interval = "right"), ll_mid)))
  expect_error(betaregscale:::brs_loglik(par, y ~ x, raw, ncuts = 10, interval = "zz"))
  expect_equal(
    betaregscale:::brs_loglik_var(c(par, 0), y ~ x | x, raw, ncuts = 10, interval = "m"),
    betaregscale:::brs_loglik_var(c(par, 0), y ~ x | x, raw, ncuts = 10), tolerance = 1e-12)
  expect_equal(betaregscale:::compute_start(y ~ x, raw, ncuts = 10, interval = "m"),
               betaregscale:::compute_start(y ~ x, raw, ncuts = 10), tolerance = 1e-12)
})

test_that("L4-21: advisory lim warnings come once from the parent call, not per refit", {
  n_adv <- function(expr) {
    k <- 0L
    withCallingHandlers(expr, warning = function(w) {
      if (grepl("cover only|is ignored for interval", conditionMessage(w))) k <<- k + 1L
      invokeRestart("muffleWarning")
    })
    k
  }
  raw <- .l4_raw(.l4_sim("mid", n = 60L))
  expect_identical(n_adv(fit <- brs(y ~ x, data = raw, ncuts = 10, lim = 0.3)), 1L)
  set.seed(9L)
  expect_identical(n_adv(brs_bootstrap(fit, R = 10L)), 0L)
  set.seed(9L)
  expect_identical(n_adv(brs_bootstrap(fit, R = 10L, ci_type = "bca")), 0L)
  set.seed(9L)
  expect_identical(n_adv(brs_cv(y ~ x, data = raw, k = 3L, ncuts = 10, lim = 0.3)), 1L)
})

test_that("L4-22: cells squeezed by the 1e-5 clamp are reported, not taken for D2", {
  expect_error(brs_check(c(1, 5), ncuts = 2e5), "`ncuts` too large")
  expect_error(suppressMessages(brs_prep(data.frame(y = c(1, 5)), ncuts = 2e5)),
               "`ncuts` too large")
  expect_error(suppressMessages(brs_prep(data.frame(left = 1e-6, right = 5e-6),
                                         ncuts = 1)),
               "inside the 1e-5 border clamp")
  expect_error(suppressMessages(brs_prep(data.frame(y = 50, left = 50, right = 50),
                                         ncuts = 100)),
               "delta = 3 with left >= right")
})

test_that("L4-23: predict(type = 'score') support and its offset from the expected score", {
  s <- .l4_sim("right", n = 400L, seed = 23L)
  raw <- .l4_raw(s)
  fr <- brs(y ~ x, data = raw, ncuts = 10, interval = "right")
  fl <- brs(y ~ x, data = raw, ncuts = 10, interval = "left")
  sr <- predict(fr, type = "score")
  sl <- predict(fl, type = "score")
  es <- predict(fr, type = "expected_score")
  expect_true(all(sr > 0 & sr < 11) && all(sl > -1 & sl < 10))
  expect_lt(abs(mean(sr - es) - 0.5), 0.05)
  expect_lt(abs(mean(sl - es) + 0.5), 0.05)
})

test_that("L4-24: a row with only delta (no score, no bounds) is rejected for every delta", {
  for (dl in 0:3) {
    d <- data.frame(y = c(3, NA), left = NA_real_, right = NA_real_, delta = c(NA, dl),
                    x = 1:2)
    expect_error(suppressMessages(brs_prep(d, ncuts = 10)),
                 "Observation\\(s\\) 2: all relevant columns are NA", label = dl)
  }
  # the zero-width check tolerates NA endpoints instead of failing cryptically
  expect_silent(betaregscale:::.brs_stop_zero_width(c(FALSE, NA), c(1, NA), c(TRUE, TRUE)))
})
