# ============================================================================ #
# Lote 7 (2026-09 audit): bootstrap DGP, post-fit diagnostics, brsmm inference,
# low findings. Named cases L7-*. Recovery checks use 3-4 Monte Carlo SE.
# ============================================================================ #

.l7_eps <- 1e-5

# Scores on 0..K from a beta regression with a factor and a positive covariate.
.l7_data <- function(n = 120L, seed = 71L, K = 20L) {
  set.seed(seed)
  x <- runif(n, 0.5, 3)
  g <- factor(sample(c("a", "b", "c"), n, TRUE))
  z <- rnorm(n)
  mu <- plogis(-0.3 + 0.6 * log(x) + c(a = 0, b = 0.4, c = -0.4)[as.character(g)])
  sh <- brs_repar(mu, 0.2, 2L)
  ys <- rbeta(n, sh$shape1, sh$shape2)
  data.frame(y = round(ys * K), pain = round(ys * K), x = x, g = g, z = z,
             yu = pmin(pmax(ys, 1e-4), 1 - 1e-4))
}

# se_boot / Wald SE for a fit (R replicates, fixed seed).
.l7_ratio <- function(fit, R = 40L) {
  set.seed(5)
  b <- brs_bootstrap(fit, R = R)
  list(boot = b, ratio = b$se_boot / sqrt(diag(vcov(fit))))
}

# L7-A: parametric bootstrap --------------------------------------------------

test_that("L7-A3a: bootstrap refits keep the fit's formula (any response name, factor, log, no intercept)", {
  d <- .l7_data()
  fits <- list(
    pain = brs(pain ~ x, data = d, ncuts = 20),
    factor = brs(y ~ g, data = d, ncuts = 20),
    log = brs(y ~ log(x), data = d, ncuts = 20),
    noint = brs(y ~ 0 + x, data = d, ncuts = 20),
    var = brs(y ~ log(x) + g | z, data = d, ncuts = 20)
  )
  for (nm in names(fits)) {
    r <- .l7_ratio(fits[[nm]])
    expect_identical(attr(r$boot, "n_success"), 40L, info = nm)
    expect_identical(r$boot$parameter, names(coef(fits[[nm]])), info = nm)
    # Before: 0 replicates (error); now the bootstrap SE tracks the Wald SE
    expect_true(all(r$ratio > 0.5 & r$ratio < 2), info = nm)
  }
})

test_that("L7-A3b: a response also present in the global environment is still resimulated", {
  d <- .l7_data()
  fit <- brs(pain ~ x, data = d, ncuts = 20)
  pain <- d$pain + 0  # the old code silently reused this vector: se_boot = 0
  assign("pain", pain, envir = globalenv())
  on.exit(rm("pain", envir = globalenv()), add = TRUE)
  r <- .l7_ratio(fit, R = 20L)
  expect_true(all(r$boot$se_boot > 0))
  # Response living only in the formula environment: materialised as a column
  d2 <- d[, c("x", "g", "z")]
  fit2 <- brs(pain ~ x, data = d2, ncuts = 20)
  set.seed(9)
  b2 <- brs_bootstrap(fit2, R = 12L)
  expect_true(all(b2$se_boot > 0))
})

test_that("L7-A3c: an expression on the left-hand side is refused with a clear message", {
  d <- .l7_data()
  fit <- brs(I(y + 0) ~ x, data = d, ncuts = 20)
  expect_error(brs_bootstrap(fit, R = 10L), "needs the response to be a variable")
})

test_that("L7-A4a: exact data stay exact (delta = 0) in every replicate", {
  d <- .l7_data()
  fit <- suppressMessages(brs(yu ~ x, data = d))
  expect_true(all(fit$delta == 0L))
  st <- betaregscale:::.brs_boot_setup(fit)
  expect_true(all(st$mech == "exact"))
  set.seed(1)
  dr <- betaregscale:::.brs_boot_data(st, fit)
  expect_true(all(dr$yu > 0 & dr$yu < 1))
  fr <- betaregscale:::.brs_refit(fit, dr)
  # Before: brs_sim re-gridded them (delta = 3 for all rows)
  expect_true(all(fr$delta == 0L))
  expect_false(isTRUE(all.equal(dr$yu, d$yu)))
})

test_that("L7-A4b: scores are re-coarsened on the fit's grid (right: floor over K + 1)", {
  d <- .l7_data()
  for (iv in c("mid", "right")) {
    fit <- brs(y ~ x, data = d, ncuts = 20, interval = iv)
    st <- betaregscale:::.brs_boot_setup(fit)
    set.seed(2)
    dr <- betaregscale:::.brs_boot_data(st, fit)
    expect_true(all(dr$y %in% 0:20), info = iv)
    set.seed(2)
    ys <- rbeta(nrow(d), st$shape1, st$shape2)
    yc <- pmin(pmax(ys, .l7_eps), 1 - .l7_eps)
    s_exp <- if (iv == "mid") round(20 * yc) else floor(21 * yc)
    expect_identical(as.numeric(dr$y), as.numeric(pmin(pmax(s_exp, 0), 20)), info = iv)
  }
})

test_that("L7-A4c: analyst thresholds (brs_prep Modes 2-3) are kept; y* chooses the side", {
  set.seed(3)
  n <- 80L
  x <- runif(n)
  lo <- sample(2:8, n, TRUE)
  raw <- data.frame(left = lo, right = lo + 3, x = x)
  p <- suppressMessages(brs_prep(raw, ncuts = 20))
  fit <- brs(y ~ x, data = p)
  st <- betaregscale:::.brs_boot_setup(fit)
  expect_true(all(st$mech == "analyst"))
  l0 <- p$left
  u0 <- p$right
  seen <- integer(0)
  for (s in 1:5) {
    set.seed(s)
    dr <- betaregscale:::.brs_boot_data(st, fit)
    ok <- (dr$delta == 1L & dr$left == .l7_eps & dr$right == l0) |
      (dr$delta == 3L & dr$left == l0 & dr$right == u0) |
      (dr$delta == 2L & dr$left == u0 & dr$right == 1 - .l7_eps)
    expect_true(all(ok))
    seen <- union(seen, unique(dr$delta))
  }
  # All three outcomes occur: the rows are not frozen copies
  expect_setequal(seen, 1:3)
  # Mode 2: delta = 1 on a score keeps its upper threshold
  p2 <- suppressMessages(suppressWarnings(
    brs_prep(data.frame(y = c(5, 12, 7, 9, 3, 15, 10, 6, 8, 11) + 0,
                        delta = c(1, NA, NA, NA, NA, NA, NA, NA, NA, NA),
                        x = seq(0, 1, length.out = 10)), ncuts = 20)
  ))
  f2 <- brs(y ~ x, data = p2)
  st2 <- betaregscale:::.brs_boot_setup(f2)
  expect_identical(st2$mech, c("analyst", rep("score", 9)))
  set.seed(4)
  dr2 <- betaregscale:::.brs_boot_data(st2, f2)
  expect_true((dr2$delta[1] == 1L && dr2$right[1] == p2$right[1]) ||
                (dr2$delta[1] == 2L && dr2$left[1] == p2$right[1]))
})

test_that("L7-A4d: mixed raw input keeps exact rows in (0, 1) and scores on the grid", {
  d <- .l7_data()
  d$ym <- ifelse(seq_len(nrow(d)) %% 3 == 0, d$yu, d$y)
  fit <- suppressWarnings(brs(ym ~ x, data = d, ncuts = 20))
  st <- betaregscale:::.brs_boot_setup(fit)
  set.seed(6)
  dr <- betaregscale:::.brs_boot_data(st, fit)
  ex <- st$mech == "exact"
  expect_true(all(dr$ym[ex] > 0 & dr$ym[ex] < 1))
  expect_true(all(dr$ym[!ex] %in% 0:20))
})

test_that("L7-A5: failed replicates are counted and printed", {
  d <- .l7_data()
  fit <- brs(y ~ x, data = d, ncuts = 20)
  set.seed(7)
  b <- brs_bootstrap(fit, R = 12L)
  expect_identical(attr(b, "n_failed"), attr(b, "n_attempted") - attr(b, "n_success"))
  expect_equal(attr(b, "fail_rate"), attr(b, "n_failed") / attr(b, "n_attempted"))
  expect_output(print(b), "Failed replicates: 0 \\(0.0% of attempts\\)")
})

test_that("L7-A6: brs_sim warns when every observation is censored on one side", {
  d <- data.frame(x = seq(0, 1, length.out = 30))
  expect_warning(brs_sim(~ x, data = d, beta = c(0, 1), phi = -1, delta = 1),
                 "informative censoring")
  expect_warning(brs_sim(~ x, data = d, beta = c(0, 1), phi = -1, delta = 2),
                 "informative censoring")
  # Natural one-sided data (all scores 0): no finite MLE either
  expect_warning(brs_sim(~ x, data = d, beta = c(-14, 0), phi = -3, ncuts = 5),
                 "censored on the same side")
  expect_silent(brs_sim(~ x, data = d, beta = c(0, 1), phi = -1))
})

test_that("L7-A7: basic-interval MCSE uses the reflected quantiles", {
  d <- .l7_data()
  fit <- brs(y ~ x, data = d, ncuts = 20)
  set.seed(8)
  b <- brs_bootstrap(fit, R = 40L, ci_type = "basic", keep_draws = TRUE)
  dr <- attr(b, "boot_draws")
  m <- betaregscale:::.boot_mcse_limits(dr[, 2], probs = c(0.025, 0.975))
  # Lower basic limit 2 theta - q(0.975): its MCSE is that of q(0.975)
  expect_equal(b$mcse_lower[2], m[2])
  expect_equal(b$mcse_upper[2], m[1])
})

test_that("L7-A8: BCa warns once per session and reports jackknife failures", {
  d <- .l7_data(n = 30L)
  fit <- brs(y ~ x, data = d, ncuts = 20)
  env <- betaregscale:::.brs_once
  if (exists("bca_approx", envir = env)) rm("bca_approx", envir = env)
  set.seed(9)
  expect_warning(b <- brs_bootstrap(fit, R = 12L, ci_type = "bca"), "approximation")
  expect_identical(attr(b, "n_jack_failed"), 0L)
  expect_output(print(b), "Failed jackknife refits: 0")
  set.seed(9)
  expect_no_warning(brs_bootstrap(fit, R = 12L, ci_type = "bca"))
})

# L7-B: post-fit diagnostics --------------------------------------------------

test_that("L7-B1: all left-censored data: clamp and Hessian warnings, SEs NA (not 0)", {
  d0 <- data.frame(y = rep(0, 30), x = seq(-1, 1, length.out = 30))
  w <- character(0)
  f0 <- withCallingHandlers(brs(y ~ 1, data = d0, ncuts = 10),
                            warning = function(w_) {
                              w <<- c(w, conditionMessage(w_))
                              invokeRestart("muffleWarning")
                            })
  expect_true(any(grepl("30 of 30 observations .* clamp boundary", w)))
  expect_true(any(grepl("Hessian near-singular or not negative definite", w)))
  expect_false(f0$diagnostics$hessian_nd)
  expect_identical(f0$diagnostics$n_clamped, 30L)
  expect_warning(s <- summary(f0), "set to NA|singular")
  expect_true(is.na(s$coefficients$mean[1, "Std. Error"]))
})

test_that("L7-B2: near-collinear and aliased columns are flagged before optim", {
  set.seed(10)
  x <- rnorm(150)
  s <- suppressWarnings(brs_sim(~ x, data = data.frame(x = x), beta = c(0.2, 0.5),
                                phi = -1.5, ncuts = 10))
  s$x3 <- x + rnorm(150, 0, 1e-5)
  s$x4 <- 2 * s$x - 1
  # The exact Hessian of this design is near-singular too: both warnings, none leaks
  w <- character(0)
  withCallingHandlers(brs(y ~ x + x3, data = s), warning = function(w_) {
    w <<- c(w, conditionMessage(w_))
    invokeRestart("muffleWarning")
  })
  expect_true(any(grepl("nearly collinear.*'x', 'x3'", w)))
  expect_true(any(grepl("Hessian near-singular", w)))
  expect_error(brs(y ~ x + x4, data = s), "rank deficient.*'x4'")
  s$z4 <- s$x4
  expect_error(brs(y ~ x | x + z4, data = s), "precision model matrix is rank deficient")
})

test_that("L7-B3: a regular fit is silent and its diagnostics are clean", {
  set.seed(11)
  d <- data.frame(x = rnorm(200))
  s <- brs_sim(~ x, data = d, beta = c(0.2, 0.5), phi = -1.5, ncuts = 10)
  expect_silent(fit <- brs(y ~ x, data = s))
  dg <- fit$diagnostics
  expect_true(dg$hessian_nd)
  expect_lt(dg$grad_step, 0.05)
  expect_identical(dg$n_clamped, 0L)
  # The stored gradient is the compiled one: it agrees with numDeriv at the estimate
  ll <- function(p) betaregscale:::.brs_loglik_fixed_cpp(
    p, fit$model_matrices$X, fit$Y[, "left"], fit$Y[, "right"], fit$Y[, "yt"],
    as.integer(fit$delta), 0L, 0L, 2L)
  g_nd <- numDeriv::grad(ll, unname(fit$par))
  expect_equal(dg$grad_norm, max(abs(g_nd)), tolerance = 1e-3)
  # And the central-difference helper (brsmm) agrees with numDeriv on a smooth function
  expect_equal(betaregscale:::.brs_num_grad(ll, unname(fit$par) + 0.1),
               numDeriv::grad(ll, unname(fit$par) + 0.1), tolerance = 1e-5)
})

test_that("L7-B4: vcov never uses a generalised inverse", {
  set.seed(12)
  d <- data.frame(x = rnorm(60))
  s <- brs_sim(~ x, data = d, beta = c(0.2, 0.5), phi = -1.5, ncuts = 10)
  fit <- brs(y ~ x, data = s)
  bad <- fit
  bad$hessian[] <- 0
  expect_warning(V <- vcov(bad), "singular")
  expect_true(all(is.na(V)))
  ind <- fit
  ind$hessian <- -diag(c(1, -1, 1))
  dimnames(ind$hessian) <- dimnames(fit$hessian)
  expect_warning(V2 <- vcov(ind), "1 negative or non-finite variance")
  expect_true(all(is.na(V2[2, ])) && all(is.na(V2[, 2])))
  expect_equal(unname(diag(V2)[c(1, 3)]), c(1, 1))
})

test_that("L7-B5: brsmm with sigma_b = 0 reports a variance component at the boundary", {
  set.seed(4)
  G <- 30L
  m <- 6L
  xx <- rnorm(G * m)
  sm <- suppressWarnings(brs_sim(~ xx, data = data.frame(xx = xx), beta = c(0.1, 0.4),
                                 phi = -1.5, ncuts = 10))
  sm$id <- factor(rep(seq_len(G), each = m))
  expect_warning(fm <- brsmm(y ~ xx, random = ~ 1 | id, data = sm),
                 "Variance component at the boundary")
  expect_true(fm$diagnostics$re_boundary)
  expect_lt(fm$diagnostics$re_gain, 1e-3)
  # Rank-deficient random-effects design
  sm$xx2 <- 2 * sm$xx
  expect_error(brsmm(y ~ xx, random = ~ 1 + xx + xx2 | id, data = sm),
               "random-effects model matrix is rank deficient")
})

# L7-C: brsmm inference -------------------------------------------------------

test_that("L7-C1: ICC of logit(Y) matches a Monte Carlo value (repar 2, 1 and 0)", {
  mc <- function(b0, sb, phi, repar, link, G = 60000L) {
    set.seed(123)
    b <- rnorm(G, 0, sb)
    mu <- betaregscale:::.clamp_mu_by_repar(betaregscale:::apply_inv_link(b0 + b, link), repar)
    sh <- brs_repar(mu, rep(phi, G), repar)
    l1 <- qlogis(rbeta(G, sh$shape1, sh$shape2))
    l2 <- qlogis(rbeta(G, sh$shape1, sh$shape2))
    r <- cor(l1, l2)
    c(r, (1 - r^2) / sqrt(G))
  }
  cases <- list(c(0.3, 0.6, 0.2, 2), c(0.5, 0.3, 30, 1), c(log(2), 0.5, 3, 0))
  for (cs in cases) {
    link <- if (cs[4] == 0) "log" else "logit"
    icc <- betaregscale:::.brs_icc_logit(cs[1], cs[3], cs[2]^2, link, cs[4])
    m <- mc(cs[1], cs[2], cs[3], cs[4], link)
    expect_lt(abs(icc - m[1]), 4 * m[2])
  }
  # The logistic formula is far off for a beta response (0.0986 vs ~0.30 here)
  expect_gt(betaregscale:::.brs_icc_logit(0.3, 0.2, 0.36, "logit", 2L) - 0.36 / (0.36 + pi^2 / 3), 0.15)
})

test_that("L7-C2: brsmm_re_study uses the beta ICC and names the terms", {
  set.seed(13)
  G <- 12L
  id <- factor(rep(seq_len(G), each = 8))
  x1 <- rnorm(length(id))
  mu <- plogis(0.2 + 0.5 * x1 + rnorm(G, sd = 0.6)[id])
  sh <- brs_repar(mu, 0.3, 2L)
  d <- data.frame(y = round(rbeta(length(id), sh$shape1, sh$shape2) * 20), x1 = x1, id = id)
  fm <- suppressWarnings(brsmm(y ~ x1, random = ~ 1 | id, data = d, ncuts = 20))
  rs <- brsmm_re_study(fm)
  eta0 <- drop(fm$model_matrices$X %*% fm$coefficients$mean)
  phi <- plogis(fm$coefficients$precision[[1]])
  expect_equal(rs$icc, betaregscale:::.brs_icc_logit(eta0, phi, fm$random$sd_b^2, "logit", 2L))
  expect_identical(colnames(rs$D), "(Intercept)")
  expect_output(print(rs), "ICC \\(logit\\(Y\\) scale")
  expect_equal(betaregscale:::.brs_icc_logit(eta0, phi, 0, "logit", 2L), 0)
})

test_that("L7-C3: summary.brsmm reports SD/Corr with transformed CIs and no tests", {
  for (q in 1:3) {
    set.seed(q)
    th <- rnorm(q * (q + 1) / 2, 0, 0.7)
    Ja <- betaregscale:::.brsmm_varcorr_jacobian(th, q)
    Jn <- numDeriv::jacobian(function(t) betaregscale:::.brsmm_varcorr_transform(t, q), th)
    expect_equal(Ja, Jn, tolerance = 1e-7, info = paste("q =", q))
  }
  set.seed(14)
  G <- 12L
  id <- factor(rep(seq_len(G), each = 8))
  x1 <- rnorm(length(id))
  mu <- plogis(0.2 + 0.5 * x1 + rnorm(G, sd = 0.6)[id])
  sh <- brs_repar(mu, 0.3, 2L)
  d <- data.frame(y = round(rbeta(length(id), sh$shape1, sh$shape2) * 20), x1 = x1, id = id)
  fm <- suppressWarnings(brsmm(y ~ x1, random = ~ 1 | id, data = d, ncuts = 20))
  s <- summary(fm)
  expect_identical(colnames(s$coefficients$random), c("Estimate", "Std. Error"))
  th <- fm$par[[4]]
  se <- sqrt(vcov(fm)[4, 4])
  expect_equal(c(s$varcorr$lower, s$varcorr$upper), exp(th + c(-1, 1) * qnorm(0.975) * se))
  expect_equal(s$varcorr$estimate, fm$random$sd_b[[1]])
  out <- capture.output(print(s))
  expect_true(any(grepl("Random effects \\(SD and Corr", out)))
  expect_false(any(grepl("logSD", out)))
})

test_that("L7-C4: anova(brs, brsmm) uses the chi-bar-square mixture and prints why", {
  set.seed(15)
  G <- 15L
  id <- factor(rep(seq_len(G), each = 8))
  x1 <- rnorm(length(id))
  mu <- plogis(0.2 + 0.5 * x1 + rnorm(G, sd = 0.7)[id])
  sh <- brs_repar(mu, 0.3, 2L)
  d <- data.frame(y = round(rbeta(length(id), sh$shape1, sh$shape2) * 20), x1 = x1, id = id)
  m0 <- brs(y ~ x1, data = d, ncuts = 20)
  m1 <- suppressWarnings(brsmm(y ~ x1, random = ~ 1 | id, data = d, ncuts = 20))
  tab <- anova(m0, m1)
  lr <- tab$Chisq[2]
  expect_equal(tab$`Pr(>Chisq)`[2], 0.5 * pchisq(lr, 1, lower.tail = FALSE))
  expect_output(print(tab), "chi-bar-square mixture")
  # Naive chi2(1) p-value is twice as large
  expect_equal(pchisq(lr, 1, lower.tail = FALSE), 2 * tab$`Pr(>Chisq)`[2])
})

test_that("L7-C5: BIC.brsmm counts observations", {
  set.seed(16)
  G <- 10L
  id <- factor(rep(seq_len(G), each = 6))
  x1 <- rnorm(length(id))
  sh <- brs_repar(plogis(0.2 + 0.5 * x1 + rnorm(G, sd = 0.5)[id]), 0.3, 2L)
  d <- data.frame(y = round(rbeta(length(id), sh$shape1, sh$shape2) * 20), x1 = x1, id = id)
  fm <- suppressWarnings(brsmm(y ~ x1, random = ~ 1 | id, data = d, ncuts = 20))
  expect_equal(BIC(fm), -2 * fm$value + log(nobs(fm)) * fm$npar)
})

# L7-D: low findings ----------------------------------------------------------

test_that("L7-D1: brsmm control entries are merged into the default maxit", {
  mc <- betaregscale:::.brs_merge_control
  expect_identical(mc(list(maxit = 2000L), list(reltol = 1e-10)),
                   list(maxit = 2000L, reltol = 1e-10))
  expect_identical(mc(list(maxit = 2000L), list(maxit = 5L)), list(maxit = 5L))
  expect_identical(mc(list(maxit = 2000L), NULL), list(maxit = 2000L))
  expect_error(mc(list(maxit = 2000L), list(5)), "named list")
})

test_that("L7-D2: summary() leaves the RNG state untouched (brs and brsmm)", {
  set.seed(17)
  d <- data.frame(x = rnorm(60), id = factor(rep(1:6, each = 10)))
  s <- brs_sim(~ x, data = d, beta = c(0.2, 0.5), phi = -1.5, ncuts = 10)
  s$id <- d$id
  fit <- brs(y ~ x, data = s)
  set.seed(1)
  invisible(summary(fit))
  r1 <- runif(3)
  set.seed(1)
  r2 <- runif(3)
  expect_identical(r1, r2)
  fm <- suppressWarnings(brsmm(y ~ x, random = ~ 1 | id, data = s))
  set.seed(2)
  invisible(summary(fm))
  r3 <- runif(3)
  set.seed(2)
  expect_identical(r3, runif(3))
})

test_that("L7-D3: marginal effects use the precision terms, not q (| 0 + z)", {
  set.seed(18)
  d <- data.frame(x = rnorm(150), z = runif(150, 0.5, 2))
  s <- brs_sim(~ x | 0 + z, data = d, beta = c(0.2, 0.5), zeta = -1, ncuts = 10)
  # brs_sim keeps z (it used to drop the first column of Z, the intercept or not)
  expect_equal(s$z, d$z)
  fit <- brs(y ~ x | 0 + z, data = s)
  expect_identical(fit$q, 1L)
  ph <- betaregscale:::.brs_me_predict(fit, s, unname(fit$par), "precision", "response")
  expect_equal(ph, fitted(fit, type = "phi"), tolerance = 1e-12)
  me <- brs_marginaleffects(fit, model = "precision", interval = FALSE)
  expect_identical(me$variable, "z")
  expect_true(is.finite(me$ame))
})

test_that("L7-D4: one per-observation rule for (0, 1) values mixed with scores", {
  expect_warning(r <- brs_check(c(0.3, 5, 10), ncuts = 10), "mixes values in \\(0, 1\\)")
  # Before: 0.3 was read as a score (delta = 3, cell [-0.02, 0.08] / 10)
  expect_identical(unname(r[, "delta"]), c(0, 3, 2))
  expect_equal(unname(r[1, c("left", "right", "yt")]), c(0.3, 0.3, 0.3))
  expect_warning(p <- suppressMessages(brs_prep(data.frame(y = c(0.3, 5, 10)), ncuts = 10)),
                 "mixes values in \\(0, 1\\)")
  expect_identical(as.numeric(p$delta), as.numeric(r[, "delta"]))
  expect_equal(p$left, unname(r[, "left"]))
  expect_equal(p$right, unname(r[, "right"]))
  # 0 is a score: (0, 1) values next to 0 only are not ambiguous
  expect_no_warning(r0 <- brs_check(c(0, 0.5), ncuts = 10))
  expect_identical(unname(r0[, "delta"]), c(1, 0))
})

test_that("L7-D5: repar 2 variable-dispersion start: moment intercept, zero slopes", {
  set.seed(19)
  d <- data.frame(x = rnorm(80), z = rnorm(80))
  s <- brs_sim(~ x | z, data = d, beta = c(2.5, 0.2), zeta = c(-2, 0.3), ncuts = 10)
  st <- compute_start(y ~ x | z, data = s, ncuts = 10)
  expect_identical(unname(st[["phi_z"]]), 0)
  st1 <- compute_start(y ~ x, data = s, ncuts = 10)
  expect_equal(unname(st[["phi_(Intercept)"]]), unname(st1[["phi"]]))
})

test_that("L7-D6: brs_prep warns on rows whose bounds cover the whole scale", {
  raw <- data.frame(left = c(2, -0.5), right = c(5, 10.5), x = 1:2)
  expect_warning(p <- suppressMessages(brs_prep(raw, ncuts = 10)),
                 "Observation\\(s\\) 2: the interval covers the whole scale")
  expect_identical(p$delta, c(3L, 3L))
})

# L7-V: validator changes (theory C1/C2, numerics C1, cheap items) -----------

test_that("L7-V1: a one-sided analyst row keeps its threshold; delta follows y*", {
  set.seed(31)
  n <- 60L
  raw <- data.frame(right = sample(3:7, n, TRUE), x = runif(n))
  p <- suppressMessages(brs_prep(raw, ncuts = 10))
  expect_true(all(p$delta == 1L))
  fit <- suppressWarnings(brs(y ~ x, data = p))
  st <- betaregscale:::.brs_boot_setup(fit)
  expect_true(all(st$mech == "analyst"))
  c0 <- p$right
  set.seed(32)
  dr <- betaregscale:::.brs_boot_data(st, fit)
  set.seed(32)
  ys <- pmin(pmax(rbeta(n, st$shape1, st$shape2), .l7_eps), 1 - .l7_eps)
  # y* <= c -> delta 1 on [eps, c]; else delta 2 on [c, 1 - eps]
  expect_identical(dr$delta, ifelse(ys <= c0, 1L, 2L))
  expect_equal(dr$right[dr$delta == 1L], c0[dr$delta == 1L])
  expect_equal(dr$left[dr$delta == 2L], c0[dr$delta == 2L])
  expect_true(all(dr$left[dr$delta == 1L] == .l7_eps & dr$right[dr$delta == 2L] == 1 - .l7_eps))
})

test_that("L7-V2: ICC is NA with a warning when the clamp of the mean drives it", {
  icc <- betaregscale:::.brs_icc_logit
  # logit: finite moments, no warning
  expect_no_warning(v <- icc(0.3, 0.2, 0.36, "logit", 2L))
  expect_true(is.finite(v))
  # probit with sigma_b^2 >= 1/2 and cloglog with a large sigma_b: infinite moments
  expect_warning(v1 <- icc(0.2, 0.1, 0.64, "probit", 2L), "ICC not available")
  expect_true(is.na(v1))
  expect_warning(v2 <- icc(0.2, 0.2, 1, "cloglog", 2L), "driven by the clamp")
  expect_true(is.na(v2))
  # small sigma_b with the same links: moments dominated by the bulk, value kept
  expect_no_warning(icc(0.2, 0.1, 0.09, "probit", 2L))
  expect_no_warning(icc(0.2, 0.2, 0.09, "cloglog", 2L))
  # through brsmm_re_study on a probit fit with a large random intercept
  set.seed(12)
  G <- 40
  id <- factor(rep(1:G, each = 8))
  x1 <- rnorm(length(id))
  sh <- brs_repar(pnorm(0.2 + 0.4 * x1 + rnorm(G, sd = 1.2)[id]), 0.2, 2L)
  d <- data.frame(y = round(rbeta(length(id), sh$shape1, sh$shape2) * 20), x1 = x1, id = id)
  fp <- suppressWarnings(brsmm(y ~ x1, random = ~ 1 | id, data = d, ncuts = 20, link = "probit"))
  expect_warning(rs <- brsmm_re_study(fp), "ICC not available")
  expect_true(is.na(rs$icc))
})

test_that("L7-V3: the gradient check is the log-likelihood gain of a Newton step", {
  # Healthy fit: gain far below 0.01 and no warning
  set.seed(11)
  s <- brs_sim(~ x, data = data.frame(x = rnorm(200)), beta = c(0.2, 0.5),
               phi = -1.5, ncuts = 10)
  expect_silent(fit <- brs(y ~ x, data = s))
  dg <- fit$diagnostics
  expect_lt(dg$grad_gain, 1e-2)
  # gain = 0.5 g' (-H)^{-1} g from the stored numbers
  g <- -betaregscale:::.brs_grad_fixed_cpp(
    unname(fit$par), fit$model_matrices$X, fit$Y[, "left"], fit$Y[, "right"],
    fit$Y[, "yt"], as.integer(fit$delta), 0L, 0L, 2L)
  expect_equal(dg$grad_gain, 0.5 * sum(g * solve(-fit$hessian, g)), tolerance = 1e-8)
  # Badly scaled cubic term (x in [1e5, 1e5 + 10]): the fit stops short of the
  # optimum and the user is told to rescale
  set.seed(1)
  x <- runif(200, 1e5, 1e5 + 10)
  s2 <- suppressWarnings(brs_sim(~ xs, data = data.frame(xs = (x - 1e5) / 10),
                                 beta = c(-0.5, 1), phi = qlogis(0.2), ncuts = 10))
  s2$x <- x
  w <- character(0)
  f3 <- withCallingHandlers(brs(y ~ I(x^3), data = s2), warning = function(w_) {
    w <<- c(w, conditionMessage(w_))
    invokeRestart("muffleWarning")
  })
  # Direct check: the same model on a standardised covariate reaches the optimum
  # (-435.216), 0.446 above the I(x^3) fit (Lote 5: the former grad_gain check here
  # passed only through a wrong numDeriv Hessian at this scale, gain 3.2e17)
  fs <- suppressWarnings(brs(y ~ scale(x^3), data = s2))
  expect_gt(as.numeric(logLik(fs) - logLik(f3)), 0.1)
  # The exact (cpp) Hessian shows the real problem: -H near-singular (intercept
  # and x^3 nearly collinear), flagged together with the collinearity and the
  # advice to rescale
  expect_false(f3$diagnostics$hessian_nd)
  expect_true(any(grepl("Hessian near-singular.*rescale the covariates", w)))
  expect_true(any(grepl("nearly collinear", w)))
})

test_that("L7-V4: brs() stores the log-likelihood exactly at the returned estimate", {
  set.seed(12)
  s <- brs_sim(~ x, data = data.frame(x = rnorm(100)), beta = c(0.2, 0.5),
               phi = -1.5, ncuts = 10)
  for (m in c("BFGS", "L-BFGS-B")) {
    f <- brs(y ~ x, data = s, method = m)
    ll <- betaregscale:::.brs_loglik_fixed_cpp(
      unname(f$par), f$model_matrices$X, f$Y[, "left"], f$Y[, "right"], f$Y[, "yt"],
      as.integer(f$delta), 0L, 0L, 2L)
    expect_identical(f$value, ll, info = m)
    expect_identical(as.numeric(logLik(f)), ll, info = m)
  }
})

test_that("L7-V5: a negative gain over sd = 0 is reported as a fit that stopped short", {
  set.seed(4)
  G <- 30L
  m <- 6L
  xx <- rnorm(G * m)
  sm <- suppressWarnings(brs_sim(~ xx, data = data.frame(xx = xx), beta = c(0.1, 0.4),
                                 phi = -1.5, ncuts = 10))
  sm$id <- factor(rep(seq_len(G), each = m))
  expect_warning(fm <- brsmm(y ~ xx, random = ~ 1 | id, data = sm),
                 "sd ~ 0 has a higher log-likelihood")
  expect_lt(fm$diagnostics$re_gain, 0)
})

test_that("L7-V6: the mixed-values warning suggests the half-point rescaling", {
  expect_warning(brs_check(c(0.5, 2, 3.5, 4), ncuts = 5),
                 "half-point scores: use y \\* 2 and ncuts \\* 2")
})
