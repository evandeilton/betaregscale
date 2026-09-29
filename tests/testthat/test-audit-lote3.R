# ============================================================================ #
# Lote 3 (2026-09 audit): reparameterizations x links, E[Y], clamps, attributes
# Named cases L3-*. Numeric identities use 1e-10; numDeriv comparisons 1e-6
# where the likelihood is smooth (the clamps are kinks: see L3-6b).
# ============================================================================ #

.l3_design <- function(seed = 31L, n = 200L) {
  set.seed(seed)
  data.frame(x = runif(n, -1, 1), z = runif(n, -1, 1))
}

# Simulate from the DGP of each repar with the default links. The repar 0
# case uses p, q < 1 (U-shaped densities); p >= 1 is covered by L3-3b/L3-7c.
.l3_truth <- function(repar) {
  switch(as.character(repar),
    "0" = list(beta = c(log(0.5), 0.25), zeta = c(log(0.6), -0.2)),
    "1" = list(beta = c(0.3, 0.8), zeta = c(log(20), 0.4)),
    "2" = list(beta = c(0.3, 0.8), zeta = c(-1.5, 0.5))
  )
}

.l3_sim <- function(repar, seed = 31L, n = 200L, fixed = FALSE) {
  d <- .l3_design(seed, n)
  tr <- .l3_truth(repar)
  if (fixed) {
    brs_sim(y ~ x, data = d, beta = tr$beta, phi = tr$zeta[1], repar = repar)
  } else {
    brs_sim(y ~ x | z, data = d, beta = tr$beta, zeta = tr$zeta, repar = repar)
  }
}

# L3-1: link compatibility matrix and defaults --------------------------------

test_that("L3-1a: NULL links resolve by repar (0 log/log, 1 logit/log, 2 logit/logit)", {
  expect_equal(unname(betaregscale:::.brs_default_links(0L)), c("log", "log"))
  expect_equal(unname(betaregscale:::.brs_default_links(1L)), c("logit", "log"))
  expect_equal(unname(betaregscale:::.brs_default_links(2L)), c("logit", "logit"))
  for (rp in 0:2) {
    s <- .l3_sim(rp, fixed = TRUE)
    f <- brs(y ~ x, data = s, repar = rp)
    def <- betaregscale:::.brs_default_links(rp)
    expect_identical(f$link, def[["link"]])
    expect_identical(f$link_phi, def[["link_phi"]])
    expect_identical(attr(s, "is_prepared"), TRUE)
  }
})

test_that("L3-1b: incompatible link/repar combinations stop with the table", {
  s1 <- .l3_sim(1L, fixed = TRUE)
  s2 <- .l3_sim(2L, fixed = TRUE)
  s0 <- .l3_sim(0L, fixed = TRUE)
  expect_error(brs(y ~ x, data = s1, repar = 1L, link_phi = "logit"),
               "not compatible with `repar = 1`")
  expect_error(brs(y ~ x, data = s2, repar = 2L, link_phi = "log"),
               "not compatible with `repar = 2`")
  expect_error(brs(y ~ x, data = s0, repar = 0L, link = "logit"),
               "not compatible with `repar = 0`")
  expect_error(brs(y ~ x, data = s2, repar = 2L, link = "log"),
               "repar = 0 \\(shapes p, q\\)")
  for (bad in c("identity", "inverse", "1/mu^2")) {
    expect_error(brs(y ~ x, data = s1, repar = 1L, link_phi = bad),
                 "not accepted for positive parameters", fixed = TRUE)
    expect_error(brs(y ~ x, data = s0, repar = 0L, link = bad),
                 "not accepted for positive parameters", fixed = TRUE)
  }
  # exact matching: "logi" used to partially match "logit"
  expect_error(brs(y ~ x, data = s2, link = "logi"), "not compatible")
  expect_error(brs_sim(y ~ x, data = .l3_design(), beta = c(0, 0), phi = 0,
                       repar = 1L, link_phi = "cloglog"), "not compatible")
  expect_error(brsmm(y ~ x, random = ~ 1 | x, data = s1, repar = 1L,
                     link_phi = "probit"), "not compatible")
  expect_error(betaregscale:::.resolve_links(link = c("logit", "probit"), repar = 2L),
               "single link name")
})

test_that("L3-1c: every allowed combination fits and stores the resolved links", {
  tab <- betaregscale:::.brs_link_table
  d <- .l3_design(n = 120L)
  for (rp in 0:2) {
    for (lk in tab[[as.character(rp)]]$mu) {
      for (lp in tab[[as.character(rp)]]$phi) {
        beta <- switch(lk, log = c(-0.6, 0.2), sqrt = c(0.8, 0.1), c(0.3, 0.6))
        zeta <- switch(lp, log = c(2.5, 0.3), sqrt = c(4, 0.5), c(-1.5, 0.4))
        if (rp == 0L) zeta <- switch(lp, log = c(-0.4, 0.2), sqrt = c(0.85, 0.1))
        s <- brs_sim(y ~ x | z, data = d, beta = beta, zeta = zeta,
                     link = lk, link_phi = lp, repar = rp)
        f <- suppressWarnings(brs(y ~ x | z, data = s, link = lk, link_phi = lp,
                                  repar = rp))
        expect_identical(c(f$link, f$link_phi), c(lk, lp))
        expect_true(all(is.finite(coef(f))))
        expect_true(all(is.finite(fitted(f))))
      }
    }
  }
})

# L3-2: brsmm validation -------------------------------------------------------

test_that("L3-2: brsmm() validates repar, ncuts, lim and links", {
  s <- .l3_sim(2L, fixed = TRUE, n = 60L)
  # raw (non-prepared) scores: explicit ncuts/lim are not overridden by attributes
  d <- data.frame(y = round(100 * s$yt), x = s$x,
                  id = factor(rep(1:6, length.out = nrow(s))))
  expect_error(brsmm(y ~ x, random = ~ 1 | id, data = d, repar = 5),
               "`repar` must be one of 0, 1, or 2")
  expect_error(brsmm(y ~ x, random = ~ 1 | id, data = d, ncuts = 1),
               "`ncuts` must be an integer >= 2")
  expect_error(brsmm(y ~ x, random = ~ 1 | id, data = d, lim = -1),
               "`lim` must be a positive finite scalar")
  expect_error(brsmm(y ~ x, random = ~ 1 | id, data = d, link_phi = "log"),
               "not compatible with `repar = 2`")
  expect_error(brsmm(y ~ x, random = ~ 1 | id, data = as.list(d)),
               "`data` must be a data.frame")
})

# L3-3: parameter recovery and E[Y] per repar ---------------------------------

test_that("L3-3a: each repar recovers its parameters within 3 SE and E[Y] = a/(a+b)", {
  for (rp in 0:2) {
    s <- .l3_sim(rp, n = 400L)
    f <- brs(y ~ x | z, data = s, repar = rp)
    tr <- .l3_truth(rp)
    truth <- c(tr$beta, tr$zeta)
    se <- sqrt(diag(vcov(f)))
    expect_true(all(abs(coef(f) - truth) < 3 * se),
                info = paste("repar", rp, ":", paste(round((coef(f) - truth) / se, 2), collapse = " ")))
    shp <- brs_repar(f$hatmu, f$hatphi, repar = rp)
    ey <- shp$shape1 / (shp$shape1 + shp$shape2)
    expect_equal(unname(fitted(f)), unname(ey), tolerance = 1e-10)
    expect_equal(unname(predict(f, type = "response")), unname(ey), tolerance = 1e-10)
    expect_lt(abs(mean(predict(f)) - mean(ey)), 1e-10)
    expect_lt(abs(mean(predict(f)) - mean(s$yt)), 0.02)
  }
})

test_that("L3-3b: repar 0 recovers a shape p >= 1 (the compiled code used to cap p at 1 - 1e-5)", {
  set.seed(3L)
  d <- .l3_design(n = 400L)
  s <- brs_sim(y ~ x, data = d, beta = c(log(3), 0.3), phi = log(2), repar = 0L)
  f <- brs(y ~ x, data = s, repar = 0L)
  truth <- c(log(3), 0.3, log(2))
  se <- sqrt(diag(vcov(f)))
  expect_true(all(is.finite(se)) && all(se > 0))
  expect_true(all(abs(coef(f) - truth) < 3 * se),
              info = paste(round((coef(f) - truth) / se, 2), collapse = " "))
  # the MLE must not be beaten by the truth (it was: -1836.7 vs -2145.9)
  ll_truth <- betaregscale:::brs_loglik(truth, y ~ x, s, repar = 0L)
  expect_gte(as.numeric(logLik(f)), ll_truth)
  expect_true(all(f$hatmu > 1))
  expect_true(all(fitted(f) > 0 & fitted(f) < 1))
})

test_that("L3-3c: audit A2 - repar 1 default link recovers a precision of 20", {
  set.seed(42L)
  d <- .l3_design(n = 300L)
  s <- brs_sim(y ~ x, data = d, beta = c(0.3, 0.8), phi = log(20), repar = 1L)
  f <- brs(y ~ x, data = s, repar = 1L)
  expect_identical(f$link_phi, "log")
  phi_hat <- unname(f$hatphi)
  expect_lt(abs(phi_hat - 20) / 20, 0.3)
  # the old default (logit on the precision) is no longer reachable
  expect_error(brs(y ~ x, data = s, repar = 1L, link_phi = "logit"), "not compatible")
})

# L3-4: repar 1/2 results unchanged; repar 0 means everywhere -----------------

test_that("L3-4a: under repar 1 and 2 the fitted mean is hatmu itself", {
  for (rp in 1:2) {
    s <- .l3_sim(rp)
    f <- brs(y ~ x | z, data = s, repar = rp)
    expect_identical(fitted(f), f$hatmu)
    expect_identical(predict(f), f$hatmu)
    expect_identical(residuals(f), f$residuals)
    expect_equal(f$residuals, as.numeric(s$yt - f$hatmu), tolerance = 1e-10)
    # weighted residuals: same numbers as the (mu, precision) formula
    prec <- if (rp == 1L) f$hatphi else (1 - f$hatphi) / f$hatphi
    mu <- f$hatmu
    ystar <- qlogis(s$yt)
    mustar <- digamma(mu * prec) - digamma((1 - mu) * prec)
    v <- trigamma(mu * prec) + trigamma((1 - mu) * prec)
    expect_equal(unname(residuals(f, type = "weighted")),
                 unname((ystar - mustar) / sqrt(prec * v)), tolerance = 1e-10)
    expect_equal(unname(residuals(f, type = "sweighted")),
                 unname((ystar - mustar) / sqrt(v)), tolerance = 1e-10)
    # deviance: saturated shapes brs_repar(y, phi) == mean y with precision a+b
    y_safe <- pmin(pmax(s$yt, 1e-7), 1 - 1e-7)
    sh <- brs_repar(mu, f$hatphi, repar = rp)
    sat <- brs_repar(y_safe, f$hatphi, repar = rp)
    dev_old <- sign(s$yt - mu) * sqrt(abs(2 * (
      dbeta(y_safe, sat$shape1, sat$shape2, log = TRUE) -
        dbeta(y_safe, sh$shape1, sh$shape2, log = TRUE))))
    expect_equal(unname(residuals(f, type = "deviance")), unname(dev_old),
                 tolerance = 1e-8)
  }
})

test_that("L3-4b: under repar 0 every user-facing mean is p/(p+q), hatmu stays p", {
  s <- .l3_sim(0L)
  f <- brs(y ~ x | z, data = s, repar = 0L)
  shp <- brs_repar(f$hatmu, f$hatphi, repar = 0L)
  expect_equal(unname(f$hatmu), unname(shp$shape1), tolerance = 1e-10)
  ey <- shp$shape1 / (shp$shape1 + shp$shape2)
  expect_equal(unname(fitted(f)), unname(ey), tolerance = 1e-10)
  expect_equal(unname(predict(f)), unname(ey), tolerance = 1e-10)
  expect_equal(unname(predict(f, newdata = f$data)), unname(ey), tolerance = 1e-10)
  expect_equal(unname(residuals(f)), unname(s$yt - ey), tolerance = 1e-10)
  expect_equal(unname(residuals(f, type = "pearson")),
               unname((s$yt - ey) / sqrt(predict(f, type = "variance"))),
               tolerance = 1e-10)
  # the AME machinery predicts E[Y] too
  me_pred <- betaregscale:::.brs_me_predict(f, f$data, par = f$par,
                                            model = "mean", type = "response")
  expect_equal(unname(me_pred), unname(ey), tolerance = 1e-10)
  # marginal effect on the mean of the numeric covariate: finite differences of E[Y]
  me <- brs_marginaleffects(f, variables = "x", interval = FALSE)
  h <- 1e-6
  dp <- f$data; dm <- f$data
  dp$x <- dp$x + h; dm$x <- dm$x - h
  fd <- mean((predict(f, newdata = dp) - predict(f, newdata = dm)) / (2 * h))
  expect_equal(me$ame, fd, tolerance = 1e-5)
  # scoreprob with newdata feeds (p, q), not E[Y], into brs_repar()
  P <- brs_predict_scoreprob(f, newdata = f$data[1:10, ])
  expect_equal(unname(rowSums(P)), rep(1, 10), tolerance = 1e-8)
  expect_equal(P, brs_predict_scoreprob(f)[1:10, ], tolerance = 1e-10)
  pp <- betaregscale:::.brs_predict_params(f, newdata = f$data)
  expect_equal(pp$mu, unname(f$hatmu), tolerance = 1e-10)
  expect_equal(pp$phi, unname(f$hatphi), tolerance = 1e-10)
})

test_that("L3-4c: weighted residuals under repar 0 equal repar 1 with the same (a, b)", {
  s <- .l3_sim(1L)
  f1 <- brs(y ~ x | z, data = s, repar = 1L)
  sh <- brs_repar(f1$hatmu, f1$hatphi, repar = 1L)
  f0 <- f1
  f0$repar <- 0L
  f0$hatmu <- sh$shape1
  f0$hatphi <- sh$shape2
  for (type in c("pearson", "deviance", "weighted", "sweighted")) {
    expect_equal(residuals(f0, type = type), residuals(f1, type = type),
                 tolerance = 1e-10, info = type)
  }
  expect_equal(fitted(f0), fitted(f1), tolerance = 1e-10)
  expect_equal(brs_predict_scoreprob(f0), brs_predict_scoreprob(f1), tolerance = 1e-10)
})

test_that("L3-4d: repar 0 also works in brsmm (fitted mean, residuals, predict)", {
  skip_on_cran()
  set.seed(8L)
  g <- 12L
  ni <- 10L
  n <- g * ni
  id <- factor(rep(seq_len(g), each = ni))
  x <- runif(n, -1, 1)
  b <- rnorm(g, sd = 0.25)
  # p = 0.35 exp(0.15 x + b_id) stays below 1; q = 0.6
  p <- 0.35 * exp(0.15 * x + b[as.integer(id)])
  s <- data.frame(y = round(100 * rbeta(n, p, 0.6)), x = x, id = id)
  f <- brsmm(y ~ x, random = ~ 1 | id, data = s, repar = 0L,
             control = list(maxit = 300L))
  expect_identical(c(f$link, f$link_phi), c("log", "log"))
  shp <- brs_repar(f$fitted_mu, f$fitted_phi, repar = 0L)
  ey <- shp$shape1 / (shp$shape1 + shp$shape2)
  expect_equal(unname(fitted(f)), unname(ey), tolerance = 1e-10)
  expect_equal(unname(predict(f)), unname(ey), tolerance = 1e-10)
  expect_equal(unname(residuals(f)), unname(f$Y[, "yt"] - ey), tolerance = 1e-10)
  expect_true(all(is.finite(residuals(f, type = "weighted"))))
  # deviance: saturated mean y with the same precision a + b, as residuals.brs
  # (the old form was dbeta at E[Y] with the fitted shapes: 60% exact zeros)
  y_safe <- pmin(pmax(f$Y[, "yt"], 1e-7), 1 - 1e-7)
  prec <- shp$shape1 + shp$shape2
  rd_ref <- sign(f$Y[, "yt"] - ey) * sqrt(abs(2 * (
    dbeta(y_safe, y_safe * prec, (1 - y_safe) * prec, log = TRUE) -
      dbeta(y_safe, shp$shape1, shp$shape2, log = TRUE))))
  rd <- residuals(f, type = "deviance")
  expect_equal(unname(rd), unname(rd_ref), tolerance = 1e-10)
  expect_false(any(rd == 0))
  set.seed(1L)
  expect_gt(cor(abs(rd), abs(residuals(f, type = "rqr"))), 0.5)
  P <- betaregscale:::.brs_score_prob_matrix(
    mu = f$fitted_mu, phi = f$fitted_phi, repar = 0L,
    ncuts = f$ncuts, lim = f$lim, scores = 0:f$ncuts)
  expect_equal(unname(rowSums(P)), rep(1, n), tolerance = 1e-7)
})

# L3-5: ncuts/lim attributes from brs_prep() and brs_sim() --------------------

test_that("L3-5a: brs() honours attr(data, 'ncuts'/'lim') from brs_prep()", {
  set.seed(5L)
  n <- 150L
  x <- rnorm(n)
  shp <- brs_repar(plogis(0.3 + 0.8 * x), 0.2, repar = 2L)
  y <- round(10 * rbeta(n, shp$shape1, shp$shape2))
  p10 <- suppressMessages(brs_prep(data.frame(y = y, x = x), ncuts = 10, lim = 0.5))
  f <- brs(y ~ x, data = p10)
  expect_identical(as.integer(f$ncuts), 10L)
  expect_equal(unname(rowSums(brs_predict_scoreprob(f))), rep(1, n), tolerance = 1e-8)
  expect_warning(f100 <- brs(y ~ x, data = p10, ncuts = 100),
                 "differs from the value used by brs_prep\\(\\) \\(10\\)")
  expect_identical(as.integer(f100$ncuts), 10L)
  expect_equal(coef(f100), coef(f), tolerance = 1e-10)
  expect_silent(f10 <- brs(y ~ x, data = p10, ncuts = 10))
  expect_identical(as.integer(f10$ncuts), 10L)
  # lim
  p25 <- suppressMessages(brs_prep(data.frame(y = y, x = x), ncuts = 10, lim = 0.25))
  expect_warning(fl <- brs(y ~ x, data = p25, lim = 0.5), "`lim = 0.5` differs")
  expect_equal(fl$lim, 0.25)
  expect_equal(brs(y ~ x, data = p25)$lim, 0.25)
  # non-prepared data: explicit value wins, default 100 / 0.5
  fr <- suppressMessages(brs(y ~ x, data = data.frame(y = y, x = x), ncuts = 10))
  expect_identical(as.integer(fr$ncuts), 10L)
  fd <- suppressMessages(brs(y ~ x, data = data.frame(y = y * 10, x = x)))
  expect_identical(as.integer(fd$ncuts), 100L)
  expect_equal(fd$lim, 0.5)
})

test_that("L3-5b: brs_sim() output carries the attributes and feeds every fitter", {
  d <- .l3_design(n = 90L)
  s <- brs_sim(y ~ x, data = d, beta = c(0.3, 0.8), phi = -1.5, ncuts = 20)
  expect_identical(attr(s, "is_prepared"), TRUE)
  expect_equal(attr(s, "ncuts"), 20)
  expect_equal(attr(s, "lim"), 0.5)
  f <- brs(y ~ x, data = s)
  expect_identical(as.integer(f$ncuts), 20L)
  expect_equal(unname(rowSums(brs_predict_scoreprob(f))), rep(1, nrow(s)), tolerance = 1e-8)
  s$id <- factor(rep(1:9, each = 10))
  fm <- brsmm(y ~ x, random = ~ 1 | id, data = s, control = list(maxit = 200L))
  expect_identical(as.integer(fm$ncuts), 20L)
  set.seed(1L)
  cv <- brs_cv(y ~ x, data = s, k = 3L)
  expect_true(all(is.finite(cv$log_score)))
  bt <- brs_bootstrap(f, R = 10L)
  expect_identical(attr(bt, "n_success"), 10L)
})

# L3-6: gradients vs numDeriv -------------------------------------------------

.l3_grad_case <- function(rp, lk, lp, d, seed = 77L) {
  set.seed(seed)
  beta <- switch(lk, log = c(-0.6, 0.2), sqrt = c(0.8, 0.1), c(0.3, 0.6))
  zeta <- switch(lp, log = c(2.5, 0.3), sqrt = c(4, 0.5), c(-1.5, 0.4))
  if (rp == 0L) zeta <- switch(lp, log = c(-0.4, 0.2), sqrt = c(0.85, 0.1))
  s <- brs_sim(y ~ x | z, data = d, beta = beta, zeta = zeta,
               link = lk, link_phi = lp, repar = rp)
  X <- cbind(1, s$x); Z <- cbind(1, s$z)
  lc <- betaregscale:::link_to_code(lk); lpc <- betaregscale:::link_to_code(lp)
  par_f <- c(beta, zeta[1]) + rnorm(3, 0, 0.05)
  par_v <- c(beta, zeta) + rnorm(4, 0, 0.05)
  fnf <- function(p) betaregscale:::.brs_loglik_fixed_cpp(
    p, X, s$left, s$right, s$yt, s$delta, lc, lpc, rp)
  fnv <- function(p) betaregscale:::.brs_loglik_variable_cpp(
    p, X, Z, s$left, s$right, s$yt, s$delta, lc, lpc, rp)
  gf <- betaregscale:::.brs_grad_fixed_cpp(par_f, X, s$left, s$right, s$yt, s$delta, lc, lpc, rp)
  gv <- betaregscale:::.brs_grad_variable_cpp(par_v, X, Z, s$left, s$right, s$yt, s$delta, lc, lpc, rp)
  gnf <- numDeriv::grad(fnf, par_f)
  gnv <- numDeriv::grad(fnv, par_v)
  c(fixed = max(abs(gf - gnf)) / max(abs(gnf)), variable = max(abs(gv - gnv)) / max(abs(gnv)))
}

test_that("L3-6a: C++ gradients match numDeriv for every allowed (repar, link, link_phi)", {
  tab <- betaregscale:::.brs_link_table
  d <- .l3_design(seed = 202L, n = 150L)
  for (rp in 0:2) {
    for (lk in tab[[as.character(rp)]]$mu) {
      for (lp in tab[[as.character(rp)]]$phi) {
        err <- .l3_grad_case(rp, lk, lp, d)
        expect_lt(err[["fixed"]], 1e-6, label = paste("fixed", rp, lk, lp))
        expect_lt(err[["variable"]], 1e-6, label = paste("variable", rp, lk, lp))
      }
    }
  }
})

test_that("L3-6b: gradients near the clamps: smooth just inside, exactly flat outside", {
  set.seed(9L)
  d <- .l3_design(seed = 9L, n = 120L)
  s <- brs_sim(y ~ x, data = d, beta = c(0.3, 0.6), phi = -1.5)
  X <- cbind(1, s$x)
  grads <- function(par, rp, lk, lp) {
    lc <- betaregscale:::link_to_code(lk); lpc <- betaregscale:::link_to_code(lp)
    fn <- function(p) betaregscale:::.brs_loglik_fixed_cpp(
      p, X, s$left, s$right, s$yt, s$delta, lc, lpc, rp)
    list(cpp = betaregscale:::.brs_grad_fixed_cpp(par, X, s$left, s$right, s$yt, s$delta, lc, lpc, rp),
         num = numDeriv::grad(fn, par))
  }
  rel <- function(g) max(abs(g$cpp - g$num)) / max(abs(g$num))
  # inside, one order of magnitude from the clamp: smooth
  expect_lt(rel(grads(c(0.3, 0.6, qlogis(2e-5)), 2L, "logit", "logit")), 1e-6)
  expect_lt(rel(grads(c(0.3, 0.6, qlogis(1 - 2e-5)), 2L, "logit", "logit")), 1e-4)
  expect_lt(rel(grads(c(0.3, 0.6, log(1e7)), 1L, "logit", "log")), 1e-6)
  expect_lt(rel(grads(c(0.3, 0.6, 1), 1L, "logit", "sqrt")), 1e-6)
  # outside the clamp the likelihood does not depend on phi: C++ gradient is 0
  expect_identical(grads(c(0.3, 0.6, qlogis(1e-6)), 2L, "logit", "logit")$cpp[3], 0)
  expect_identical(grads(c(0.3, 0.6, qlogis(1 - 1e-6)), 2L, "logit", "logit")$cpp[3], 0)
  expect_identical(grads(c(0.3, 0.6, log(1e9)), 1L, "logit", "log")$cpp[3], 0)
  expect_identical(grads(c(0.3, 0.6, -0.1), 1L, "logit", "sqrt")$cpp[3], 0)
  # partial clamping of the mean (some rows below 1e-5) keeps beta smooth
  expect_lt(rel(grads(c(-12, 0.6, -1.5), 2L, "logit", "logit")), 1e-6)
  # At the boundary itself the likelihood has a kink: central differences that
  # straddle it (C++ h = 1e-6, numDeriv Richardson) legitimately disagree there.
})

# L3-7: clamp mirrors R vs C++ ------------------------------------------------

test_that("L3-7a: .clamp_phi_by_repar / .clamp_mu_by_repar values at the boundaries", {
  cp <- betaregscale:::.clamp_phi_by_repar
  cm <- betaregscale:::.clamp_mu_by_repar
  x <- c(-1, 0, 1e-6, 1e-5, 2e-5, 0.5, 1 - 2e-5, 1 - 1e-5, 1 - 1e-6, 1, 1e8, 1e9, Inf, -Inf, NaN, NA)
  expect_equal(cp(x, 2L), c(1e-5, 1e-5, 1e-5, 1e-5, 2e-5, 0.5, 1 - 2e-5, 1 - 1e-5,
                            1 - 1e-5, 1 - 1e-5, 1 - 1e-5, 1 - 1e-5, 1 - 1e-5, 1 - 1e-5,
                            1 - 1e-5, 1 - 1e-5), tolerance = 1e-10)
  for (rp in 0:1) {
    expect_equal(cp(x, rp), c(1e-5, 1e-5, 1e-5, 1e-5, 2e-5, 0.5, 1 - 2e-5, 1 - 1e-5,
                              1 - 1e-6, 1, 1e8, 1e8, 1e8, 1e8, 1e8, 1e8), tolerance = 1e-10)
  }
  expect_equal(cm(c(0, 1e-5, 0.5, 1 - 1e-5, 1, Inf, NaN, NA), 2L),
               c(1e-5, 1e-5, 0.5, 1 - 1e-5, 1 - 1e-5, 1 - 1e-5, 1 - 1e-5, 1 - 1e-5))
  expect_equal(cm(c(0, 1e-5, 0.5, 1, 2, 1e9, Inf, NaN), 0L),
               c(1e-5, 1e-5, 0.5, 1, 2, 1e8, 1e8, 1e8))
})

.l3_clamp_fixture <- function() {
  set.seed(7L)
  y <- round(100 * rbeta(40, 2, 3))
  Yc <- brs_check(y, ncuts = 100)
  Yc[1, ] <- c(1e-5, 0.005, 1e-5, 0, 1)
  Yc[2, ] <- c(0.995, 1 - 1e-5, 1 - 1e-5, 100, 2)
  Yc[3, ] <- c(0.4, 0.4, 0.4, 40, 0)
  Yc
}

.l3_cpp_vs_mirror <- function(rp, lk, lp, mus, phis, Yc) {
  X <- matrix(1, nrow(Yc), 1)
  worst <- 0
  for (mu in mus) {
    for (phi in phis) {
      par <- c(betaregscale:::apply_link(mu, lk), betaregscale:::apply_link(phi, lp))
      ll_c <- betaregscale:::.brs_loglik_fixed_cpp(
        par, X, Yc[, "left"], Yc[, "right"], Yc[, "yt"], as.integer(Yc[, "delta"]),
        betaregscale:::link_to_code(lk), betaregscale:::link_to_code(lp), rp)
      sh <- brs_repar(betaregscale:::.clamp_mu_by_repar(mu, rp),
                      betaregscale:::.clamp_phi_by_repar(phi, rp), rp)
      ll_r <- sum(betaregscale:::.brs_obs_loglik(
        Yc[, "delta"], Yc[, "left"], Yc[, "right"], Yc[, "yt"], sh$shape1, sh$shape2))
      worst <- max(worst, abs(ll_c - ll_r) / max(1, abs(ll_c)))
    }
  }
  worst
}

test_that("L3-7b: C++ log-likelihood equals the R mirror on a grid through the clamps", {
  Yc <- .l3_clamp_fixture()
  mus <- c(1e-7, 1e-6, 1e-5, 2e-5, 0.5, 1 - 2e-5, 1 - 1e-5, 1 - 1e-6)
  expect_lt(.l3_cpp_vs_mirror(2L, "logit", "logit", mus,
    c(1e-7, 1e-6, 1e-5, 2e-5, 0.3, 1 - 2e-5, 1 - 1e-5, 1 - 1e-6), Yc), 1e-10)
  expect_lt(.l3_cpp_vs_mirror(1L, "logit", "log", mus,
    c(1e-6, 1e-5, 2e-5, 20, 1e8 - 1, 1e8, 1e9), Yc), 1e-10)
  # repar 0 with p < 1 (p >= 1 up to the 1e8 clamp: L3-7c)
  expect_lt(.l3_cpp_vs_mirror(0L, "log", "log",
    c(1e-6, 1e-5, 2e-5, 0.5, 0.99999), c(1e-6, 1e-5, 0.5, 5, 1e8, 1e9), Yc), 1e-10)
})

test_that("L3-7c: repar 0 first-parameter clamp [1e-5, 1e8] agrees R vs C++ for p >= 1", {
  Yc <- .l3_clamp_fixture()
  expect_lt(.l3_cpp_vs_mirror(0L, "log", "log",
    c(1e-6, 1e-5, 2e-5, 0.5, 0.99999, 1, 2, 1e8, 1e9),
    c(1e-6, 1e-5, 0.5, 5, 1e8, 1e9), Yc), 1e-10)
})

# L3-8: finiteness at the clamps ---------------------------------------------

.l3_check_finite <- function(f, label) {
  expect_true(all(is.finite(fitted(f))), label = paste(label, "fitted"))
  expect_true(all(is.finite(predict(f, type = "variance"))), label = paste(label, "variance"))
  for (type in c("response", "pearson", "deviance", "weighted", "sweighted", "rqr")) {
    r <- residuals(f, type = type)
    expect_true(all(is.finite(r)), label = paste(label, type))
  }
  P <- brs_predict_scoreprob(f)
  expect_true(all(is.finite(P)), label = paste(label, "scoreprob"))
  # interior bins are floored at 1e-10 each (99 bins): sums are 1 + O(1e-8)
  expect_equal(unname(rowSums(P)), rep(1, nrow(P)), tolerance = 1e-7)
}

test_that("L3-8a: E[Y], residuals, pseudo-R2 and scoreprob stay finite at every clamp", {
  s <- .l3_sim(2L, fixed = TRUE, n = 60L)
  f <- brs(y ~ x, data = s)
  n <- nrow(s)
  cases <- list(
    list(repar = 2L, mu = 1e-5, phi = 1e-5), list(repar = 2L, mu = 1 - 1e-5, phi = 1 - 1e-5),
    list(repar = 2L, mu = 1e-5, phi = 1 - 1e-5), list(repar = 2L, mu = 0.5, phi = 1e-5),
    list(repar = 1L, mu = 1e-5, phi = 1e-5), list(repar = 1L, mu = 1 - 1e-5, phi = 1e8),
    list(repar = 1L, mu = 0.5, phi = 1e8),
    list(repar = 0L, mu = 1e-5, phi = 1e-5), list(repar = 0L, mu = 1e-5, phi = 1e8),
    list(repar = 0L, mu = 1e8, phi = 1e-5), list(repar = 0L, mu = 1e8, phi = 1e8)
  )
  for (cs in cases) {
    g <- f
    g$repar <- cs$repar
    g$hatmu <- rep(cs$mu, n)
    g$hatphi <- cs$phi
    g$residuals <- s$yt - betaregscale:::.brs_mean(g$hatmu, g$hatphi, cs$repar)
    set.seed(1L)
    .l3_check_finite(g, sprintf("repar %d mu=%g phi=%g", cs$repar, cs$mu, cs$phi))
    ey <- fitted(g)
    r2 <- betaregscale:::.brs_pseudo_r2(rep(0, n), ey, s$yt, "logit", cs$repar)
    expect_true(is.na(r2) || is.finite(r2))
  }
  # smallest shape reachable through the clamps is 1e-10: digamma/trigamma finite
  expect_true(is.finite(digamma(1e-10)) && is.finite(trigamma(1e-10)))
  expect_true(is.finite(sqrt(2e-10 * (trigamma(1e-10) + trigamma(1e-10)))))
})

test_that("L3-8b: fits with every observation censored at one border stay finite", {
  set.seed(2L)
  for (yb in c(0, 100)) {
    d <- data.frame(y = rep(yb, 30), x = rnorm(30))
    f <- suppressWarnings(suppressMessages(brs(y ~ x, data = d)))
    expect_true(all(is.finite(coef(f))))
    expect_true(all(f$hatmu >= 1e-5 & f$hatmu <= 1 - 1e-5))
    .l3_check_finite(f, paste("border", yb))
    expect_true(is.na(f$pseudo.r.squared) || is.finite(f$pseudo.r.squared))
  }
})

# L3-9: small n, intercept-only, factors, newdata levels ----------------------

test_that("L3-9a: n = 8, intercept-only and single-covariate fits", {
  set.seed(11L)
  ds <- brs_sim(y ~ x, data = data.frame(x = rnorm(8)), beta = c(0.2, 0.3), phi = -1)
  f1 <- brs(y ~ 1, data = ds)
  expect_true(is.na(f1$pseudo.r.squared))
  expect_true(all(is.finite(predict(f1, newdata = ds[1:2, ]))))
  expect_equal(length(unique(fitted(f1))), 1L)
  f2 <- brs(y ~ x, data = ds)
  expect_true(is.finite(f2$pseudo.r.squared))
  .l3_check_finite(f2, "n = 8")
  for (rp in 0:1) {
    s <- .l3_sim(rp, fixed = TRUE, n = 8L)
    st <- betaregscale:::compute_start(y ~ x, s, repar = rp)
    expect_true(all(is.finite(st)))
    expect_s3_class(brs(y ~ 1, data = s, repar = rp), "brs")
  }
})

test_that("L3-9b: factor covariates; newdata with a subset of levels or a new level", {
  set.seed(12L)
  n <- 150L
  dg <- data.frame(x = rnorm(n), g = factor(sample(c("a", "b", "c"), n, TRUE)))
  sg <- brs_sim(y ~ x + g, data = dg, beta = c(0.2, 0.3, 0.4, -0.3), phi = -1.5)
  sg$g <- dg$g
  f <- brs(y ~ x + g | g, data = sg)
  expect_identical(f$xlevels$mean$g, levels(dg$g))
  nd_factor <- data.frame(x = c(0, 1), g = factor(c("b", "b"), levels = levels(dg$g)))
  nd_chr <- data.frame(x = c(0, 1), g = c("b", "b"))
  p_ref <- plogis(coef(f)[["(Intercept)"]] + coef(f)[["x"]] * c(0, 1) + coef(f)[["gb"]])
  expect_equal(unname(predict(f, newdata = nd_factor)), unname(p_ref), tolerance = 1e-10)
  expect_equal(unname(predict(f, newdata = nd_chr)), unname(p_ref), tolerance = 1e-10)
  expect_equal(unname(rowSums(brs_predict_scoreprob(f, newdata = nd_chr))), c(1, 1),
               tolerance = 1e-8)
  pp <- betaregscale:::.brs_predict_params(f, newdata = nd_chr)
  expect_equal(pp$mu, unname(p_ref), tolerance = 1e-10)
  expect_length(pp$phi, 2L)
  expect_error(predict(f, newdata = data.frame(x = 0, g = "zz")), "new level")
  expect_error(brs_predict_scoreprob(f, newdata = data.frame(x = 0, g = "zz")), "new level")
  me <- brs_marginaleffects(f, variables = "x", interval = FALSE)
  expect_true(is.finite(me$ame))
  expect_error(betaregscale:::.brs_me_predict(f, data.frame(x = 0, g = "zz"), par = f$par,
                                             model = "mean", type = "response"), "new level")
  # brsmm keeps the levels as well
  sg$id <- factor(rep(1:15, each = 10))
  fm <- brsmm(y ~ x + g, random = ~ 1 | id, data = sg, control = list(maxit = 150L))
  expect_identical(fm$xlevels$mean$g, levels(dg$g))
  expect_error(predict(fm, newdata = data.frame(x = 0, g = "zz", id = "1")), "new level")
  expect_length(predict(fm, newdata = data.frame(x = c(0, 1), g = "b", id = "1")), 2L)
})

# L3-10: sqrt plateau, U-shaped repar 0, moment start guard -------------------

test_that("L3-10a: sqrt link warns when a fitted linear predictor is on the plateau", {
  expect_warning(betaregscale:::.warn_sqrt_plateau(c(1, 0, -1), "sqrt", "link_phi"),
                 "<= 0 for 2 observation")
  expect_silent(betaregscale:::.warn_sqrt_plateau(c(1, 0.5), "sqrt", "link_phi"))
  expect_silent(betaregscale:::.warn_sqrt_plateau(c(-1, -2), "log", "link_phi"))
  set.seed(5L)
  n <- 200L
  z <- runif(n, -1.5, 1.5)
  x <- runif(n, -1, 1)
  d <- data.frame(x = x, z = z)
  # precision from exp(-4) to exp(5): the sqrt-link line crosses zero
  s <- brs_sim(y ~ x | z, data = d, beta = c(0.2, 0.5), zeta = c(0.5, 3), repar = 1L)
  expect_warning(f <- brs(y ~ x | z, data = s, repar = 1L, link_phi = "sqrt"),
                 "link_phi = \"sqrt\"")
  eta_phi <- coef(f)[["(phi)_(Intercept)"]] + coef(f)[["(phi)_z"]] * s$z
  expect_true(any(eta_phi <= 0))
  expect_true(all(is.finite(fitted(f))))
})

test_that("L3-10b: U-shaped repar 0 (p, q < 1): start, fit, means and scoreprob", {
  set.seed(6L)
  n <- 300L
  su <- brs_sim(y ~ x, data = data.frame(x = rnorm(n)), beta = c(log(0.5), 0),
                phi = log(0.5), repar = 0L)
  st <- betaregscale:::compute_start(y ~ x, su, repar = 0L)
  expect_true(all(is.finite(st)))
  expect_lt(exp(st[[1]]), 1)
  f <- brs(y ~ x, data = su, repar = 0L)
  se <- sqrt(diag(vcov(f)))
  expect_true(all(abs(coef(f) - c(log(0.5), 0, log(0.5))) < 3 * se))
  .l3_check_finite(f, "U-shaped")
  expect_true(all(fitted(f) > 0 & fitted(f) < 1))
  # moment start guard: degenerate midpoints (f <= 0 or v = 0) still give a start
  dz <- data.frame(y = rep(0, 20), x = rnorm(20))
  st0 <- suppressMessages(betaregscale:::compute_start(y ~ x, dz, repar = 0L))
  expect_true(all(is.finite(st0)))
})

test_that("L3-10c: repar 1 starting values with the log precision link are moment-based", {
  s <- .l3_sim(1L)
  st <- betaregscale:::compute_start(y ~ x | z, s, repar = 1L)
  expect_true(all(is.finite(st)))
  expect_equal(unname(st[4]), 0)          # zero slope on the precision
  expect_gt(exp(st[3]), 1)                # precision start on (0, Inf)
  # repar 2 keeps the GLM-based start
  s2 <- .l3_sim(2L)
  st2 <- betaregscale:::compute_start(y ~ x | z, s2, repar = 2L)
  expect_true(all(is.finite(st2)))
})
