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

# C2 -------------------------------------------------------------------------- #

# Reference: interval log-probability in log-space, tail chosen by the median.
.audit_log_int <- function(l, u, mu, prec) {
  a <- mu * prec; b <- (1 - mu) * prec
  lF <- function(x, lt) pbeta(x, a, b, lower.tail = lt, log.p = TRUE)
  m <- qbeta(0.5, a, b)
  if (l >= m) { x1 <- lF(l, FALSE); x2 <- lF(u, FALSE) }
  else if (u <= m) { x1 <- lF(u, TRUE); x2 <- lF(l, TRUE) }
  else return(log(pbeta(u, a, b) - pbeta(l, a, b)))
  x1 + log1p(-exp(x2 - x1))
}
.audit_pkg_int <- function(l, u, mu, prec) {
  # repar = 1 (precision), link_mu logit (0), link_phi log (4)
  betaregscale:::.brs_loglik_fixed_cpp(
    c(qlogis(mu), log(prec)), matrix(1), l, u, (l + u) / 2, 3L, 0L, 4L, 1L
  )
}
# Count R warnings raised while evaluating expr (they are muffled).
.audit_n_warnings <- function(expr) {
  n <- 0L
  withCallingHandlers(expr, warning = function(w) {
    n <<- n + 1L
    invokeRestart("muffleWarning")
  })
  n
}
# Random mixed-censoring data for the compiled likelihood (delta 0..3).
.audit_ll_data <- function(seed, n, mu_fun, prec, extra = NULL) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n); z1 <- rnorm(n)
  mu <- mu_fun(x1, x2); a <- mu * prec; b <- (1 - mu) * prec
  y <- round(100 * rbeta(n, a, b))
  Y <- brs_check(y, ncuts = 100L)
  delta <- as.integer(Y[, "delta"])
  ex <- seq(1, n, by = 10)   # ~10% exact observations
  delta[ex] <- 0L
  Y[ex, "yt"] <- pmin(pmax(rbeta(length(ex), a[ex], b[ex]), 1e-5), 1 - 1e-5)
  delta[seq(2, n, by = 10)] <- 1L   # ~10% left-censored (Y <= right)
  delta[seq(3, n, by = 10)] <- 2L   # ~10% right-censored (Y >= left)
  d <- list(X = cbind(1, x1, x2), Z = cbind(1, z1), l = unname(Y[, "left"]),
            r = unname(Y[, "right"]), yt = unname(Y[, "yt"]), delta = delta)
  for (e in extra) {   # extra interval-censored rows at x = 0
    d$X <- rbind(d$X, c(1, 0, 0)); d$Z <- rbind(d$Z, c(1, 0))
    d$l <- c(d$l, e[1]); d$r <- c(d$r, e[2]); d$yt <- c(d$yt, mean(e))
    d$delta <- c(d$delta, 3L)
  }
  d
}

test_that("C2: tail interval probabilities are exact, not floored at 1e-15", {
  # Named cases from the audit: the old code returned log(1e-15) = -34.54
  cases <- rbind(
    c(0.495, 0.505, 0.90, 500),   # exact ~ -184.36
    c(0.120, 0.130, 0.01, 400),   # exact ~ -40.77
    c(0.140, 0.150, 0.01, 400),   # exact ~ -49.42
    c(0.085, 0.095, 0.01, 500)    # cancellation zone, exact ~ -32.09
  )
  for (k in seq_len(nrow(cases))) {
    r <- cases[k, ]
    expect_equal(.audit_pkg_int(r[1], r[2], r[3], r[4]),
                 .audit_log_int(r[1], r[2], r[3], r[4]),
                 tolerance = 1e-8)
  }
})

test_that("C2: gradient is non-zero and stable in the former cancellation zone", {
  f <- function(b0) betaregscale:::.brs_loglik_fixed_cpp(
    c(b0, log(500)), matrix(1), 0.085, 0.095, 0.09, 3L, 0L, 4L, 1L)
  fa <- function(b0) .audit_log_int(0.085, 0.095, plogis(b0), 500)
  b0 <- qlogis(0.01)
  for (h in c(1e-4, 1e-6)) {
    d_pkg <- (f(b0 + h) - f(b0 - h)) / (2 * h)
    d_ref <- (fa(b0 + h) - fa(b0 - h)) / (2 * h)
    expect_equal(d_pkg, d_ref, tolerance = 1e-4)
  }
  g <- betaregscale:::.brs_grad_fixed_cpp(
    c(b0, log(500)), matrix(1), 0.085, 0.095, 0.09, 3L, 0L, 4L, 1L)
  expect_true(abs(g[1]) > 1)
})

test_that("C2: outliers are not trimmed from the likelihood (precision fit)", {
  set.seed(42L); n <- 200L
  d <- data.frame(x1 = rnorm(n))
  mu <- plogis(-1.4 + 0.3 * d$x1); prec <- 300
  d$y <- round(100 * rbeta(n, mu * prec, (1 - mu) * prec))
  d$y[1:4] <- c(60, 70, 80, 90)
  prep <- suppressMessages(suppressWarnings(brs_prep(d, ncuts = 100L)))
  fit <- suppressWarnings(brs(y ~ x1, data = prep, repar = 1L, link_phi = "log"))
  X <- model.matrix(~ x1, prep)
  llex <- function(p) {
    m <- plogis(X %*% p[1:2]); ph <- exp(p[3])
    sum(vapply(seq_len(n), function(i)
      .audit_log_int(prep$left[i], prep$right[i], m[i], ph), 1))
  }
  expect_equal(fit$value, llex(fit$par), tolerance = 1e-6)
  ex <- optim(fit$par, function(p) -llex(p), method = "BFGS",
              control = list(reltol = 1e-12, maxit = 1000L))
  expect_equal(exp(fit$par[3]), exp(ex$par[3]), tolerance = 0.02)
  expect_true(exp(fit$par[3]) < 100)   # old floored fit gave ~328
})

test_that("C2: R mirror .brs_obs_loglik equals the compiled likelihood", {
  d <- .audit_ll_data(21L, 120L, function(x1, x2) plogis(0.3 + 0.6 * x1 - 0.4 * x2), 25)
  expect_setequal(unique(d$delta), 0:3)
  beta <- c(0.1, 0.4, -0.3)
  # (link_mu, link_phi, repar, gamma): both link families, all reparameterisations
  combos <- list(
    list("logit",   "log",      1L, c(log(20), 0.2)),
    list("probit",  "sqrt",     1L, c(4.0, 0.3)),
    list("cauchit", "identity", 0L, c(3.0, 0.5)),
    list("probit",  "log",      0L, c(log(2), -0.2)),
    list("cloglog", "logit",    2L, c(qlogis(0.15), 0.3)),
    list("logit",   "logit",    2L, c(qlogis(0.30), -0.4))
  )
  clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)
  for (cb in combos) {
    lm <- link_to_code(cb[[1]]); lp <- link_to_code(cb[[2]]); rp <- cb[[3]]
    gam <- cb[[4]]
    mu <- clamp(apply_inv_link(drop(d$X %*% beta), cb[[1]]), 1e-5, 1 - 1e-5)
    phi_hi <- if (rp == 2L) 1 - 1e-5 else 1e8
    # fixed dispersion
    phi <- clamp(apply_inv_link(gam[1], cb[[2]]), 1e-5, phi_hi)
    sh <- brs_repar(mu, phi, repar = rp)
    ll_r <- sum(betaregscale:::.brs_obs_loglik(d$delta, d$l, d$r, d$yt,
                                                sh$shape1, sh$shape2))
    ll_c <- betaregscale:::.brs_loglik_fixed_cpp(
      c(beta, gam[1]), d$X, d$l, d$r, d$yt, d$delta, lm, lp, rp)
    expect_equal(ll_c, ll_r, tolerance = 1e-10)
    # variable dispersion
    phi <- clamp(apply_inv_link(drop(d$Z %*% gam), cb[[2]]), 1e-5, phi_hi)
    sh <- brs_repar(mu, phi, repar = rp)
    ll_r <- sum(betaregscale:::.brs_obs_loglik(d$delta, d$l, d$r, d$yt,
                                                sh$shape1, sh$shape2))
    ll_c <- betaregscale:::.brs_loglik_variable_cpp(
      c(beta, gam), d$X, d$Z, d$l, d$r, d$yt, d$delta, lm, lp, rp)
    expect_equal(ll_c, ll_r, tolerance = 1e-10)
  }
  # Same non-finite policy as C++: zero-width interval and NaN shape -> -1e6
  expect_equal(betaregscale:::.brs_obs_loglik(3L, 0.3, 0.3, 0.3, 2, 5), -1e6)
  expect_equal(betaregscale:::.brs_obs_loglik(1L, 0.3, 0.3, 0.3, NaN, 5), -1e6)
  expect_equal(betaregscale:::.brs_loglik_fixed_cpp(
    c(0.4, 7), matrix(1), 0.3, 0.3, 0.3, 3L, 8L, 8L, 1L), -1e6)
})

test_that("C2: compiled gradients match numDeriv at random and tail points", {
  chk <- function(d, pf, pv, lm = 0L, lp = 4L, rp = 1L) {
    fnf <- function(p) betaregscale:::.brs_loglik_fixed_cpp(
      p, d$X, d$l, d$r, d$yt, d$delta, lm, lp, rp)
    fnv <- function(p) betaregscale:::.brs_loglik_variable_cpp(
      p, d$X, d$Z, d$l, d$r, d$yt, d$delta, lm, lp, rp)
    gf <- betaregscale:::.brs_grad_fixed_cpp(pf, d$X, d$l, d$r, d$yt, d$delta, lm, lp, rp)
    gv <- betaregscale:::.brs_grad_variable_cpp(pv, d$X, d$Z, d$l, d$r, d$yt, d$delta, lm, lp, rp)
    nf <- numDeriv::grad(fnf, pf); nv <- numDeriv::grad(fnv, pv)
    expect_true(all(abs(nf) > 1e-3) && all(abs(nv) > 1e-3))  # away from the optimum
    expect_lt(max(abs(gf - nf) / abs(nf)), 1e-5)             # componentwise relative
    expect_lt(max(abs(gv - nv) / abs(nv)), 1e-5)
    for (H in list(numDeriv::hessian(fnf, pf), numDeriv::hessian(fnv, pv))) {
      expect_true(all(is.finite(H)))
      expect_lt(max(abs(H - t(H))) / max(abs(H)), 1e-8)
    }
  }
  # three random points on ordinary data
  d1 <- .audit_ll_data(1L, 150L, function(x1, x2) plogis(0.2 + 0.5 * x1 - 0.3 * x2), 30)
  set.seed(11L)
  for (k in 1:3) {
    pf <- c(rnorm(3, 0, 0.5), log(runif(1, 5, 80)))
    pv <- c(pf[1:3], log(runif(1, 5, 80)), rnorm(1, 0, 0.3))
    chk(d1, pf, pv)
  }
  # former cancellation zone: mu ~ 0.01, precision ~ 500, intervals of width 0.01
  d2 <- .audit_ll_data(2L, 150L, function(x1, x2) plogis(qlogis(0.01) + 0.1 * x1), 500)
  p2f <- c(qlogis(0.01), 0.1, 0.05, log(500)); p2v <- c(p2f, 0.1)
  chk(d2, p2f, p2v)
  # same, plus two outliers far in the upper tail: [0.5, 0.51] (log p ~ -324,
  # plain-scale pbeta) and [0.8, 0.81] (log p ~ -776, Laplace regime)
  d3 <- .audit_ll_data(2L, 150L, function(x1, x2) plogis(qlogis(0.01) + 0.1 * x1), 500,
                       extra = list(c(0.5, 0.51), c(0.8, 0.81)))
  tails <- betaregscale:::.brs_obs_loglik(c(3L, 3L), c(0.5, 0.8), c(0.51, 0.81),
                                          c(0.505, 0.805), 5, 495)
  expect_true(tails[1] < -300 && tails[1] > -350)
  expect_true(tails[2] < -700 && tails[2] > -800)
  chk(d3, p2f, p2v)
  # other link families / reparameterisations
  chk(d1, c(0.2, 0.5, -0.3, 0.5), c(0.2, 0.5, -0.3, 0.5, 0.1), lm = 1L, lp = 5L, rp = 1L)
  chk(d1, c(0.2, 0.5, -0.3, -1.2), c(0.2, 0.5, -0.3, -1.2, 0.2), lm = 3L, lp = 0L, rp = 2L)
})

test_that("C2: no R warnings from the likelihood at extreme shapes", {
  # identity links (code 8), repar 1, X = 1: param = c(mu, phi) gives
  # a = mu * phi, b = (1 - mu) * phi, so any (a, b) with a + b <= 1e8 is reachable
  cpp1 <- function(d, l, r, yt, a, b) betaregscale:::.brs_loglik_fixed_cpp(
    c(a / (a + b), a + b), matrix(1), l, r, yt, as.integer(d), 8L, 8L, 1L)
  shapes <- c(0.01, 0.1, 0.5, 1, 2, 10, 100, 1000, 2082, 1e4, 1e5, 1e6, 1e7)
  xg <- c(1e-5, 1e-4, 0.001, 0.01, 0.05, 0.1, 0.3, 0.5, 0.675, 0.9, 0.95, 0.99, 0.999)
  g <- expand.grid(a = shapes, b = shapes, x = xg)
  m <- g$a / (g$a + g$b)
  g <- g[g$a + g$b <= 1e8 & m > 2e-5 & m < 1 - 2e-5, ]
  g$r <- pmin(g$x + 0.01, 1 - 1e-5)
  # cases that made log-scale pbeta() warn ("bpser underflow to -Inf")
  g <- rbind(g, data.frame(a = c(2082, 1e5), b = c(40, 23), x = c(0.675, 0.95),
                           r = c(0.685, 0.96)))
  yt <- (g$x + g$r) / 2
  n_penalty <- 0L
  nw <- .audit_n_warnings({
    for (d in 0:3) {
      vc <- mapply(function(a, b, x, r, yt) cpp1(d, x, r, yt, a, b),
                   g$a, g$b, g$x, g$r, yt)
      vr <- betaregscale:::.brs_obs_loglik(rep(d, nrow(g)), g$x, g$r, yt, g$a, g$b)
      expect_true(all(is.finite(vc)))
      expect_lt(max(abs(vc - vr) / pmax(1, abs(vr))), 1e-10)
      # LOG_PENALTY is exactly -1e6; genuine contributions can be far below
      # it at extreme shapes (a = 1e7 at x = 1e-5 gives ~ -1.15e8)
      n_penalty <- n_penalty + sum(vc == -1e6)
    }
    # gradient path too (2 * npar extra evaluations per point)
    for (i in seq(1, nrow(g), by = 7)) {
      gr <- betaregscale:::.brs_grad_fixed_cpp(
        c(g$a[i] / (g$a[i] + g$b[i]), g$a[i] + g$b[i]), matrix(1),
        g$x[i], g$r[i], yt[i], 3L, 8L, 8L, 1L)
      expect_true(all(is.finite(gr)))
    }
  })
  expect_equal(nw, 0L)
  expect_equal(n_penalty, 0L)   # the Laplace fallback always applied
})

test_that("C2: the Laplace fallback joins the exact tail continuously", {
  cpp1 <- function(d, l, r, yt, a, b) betaregscale:::.brs_loglik_fixed_cpp(
    c(a / (a + b), a + b), matrix(1), l, r, yt, as.integer(d), 8L, 8L, 1L)
  # upper tail (right-censored and interval) of Beta(5, 495): S(x) crosses
  # 1e-240 near x = 0.6856; below that the Laplace approximation is used
  xs <- seq(0.63, 0.74, by = 1e-4)
  S <- pbeta(xs, 5, 495, lower.tail = FALSE)
  i <- which(diff(S >= 1e-240) != 0)
  expect_length(i, 1L)
  for (d in c(2L, 3L)) {
    v <- vapply(xs, function(x) cpp1(d, x, min(x + 0.01, 1), x + 0.005, 5, 495), 1)
    expect_true(all(is.finite(v)))
    expect_true(all(diff(v) < 0))                       # strictly decreasing
    d2 <- abs(diff(v, differences = 2))
    expect_lt(max(d2[(i - 3):(i + 3)]), 10 * median(d2))  # no jump at the switch
  }
  # lower tail mirror: Beta(495, 5), left-censored at x (crossing at 0.3144)
  xs <- rev(1 - xs)
  Fx <- pbeta(xs, 495, 5)
  i <- which(diff(Fx >= 1e-240) != 0)
  expect_length(i, 1L)
  v <- vapply(xs, function(x) cpp1(1L, 1e-5, x, x, 495, 5), 1)
  expect_true(all(diff(v) > 0))
  d2 <- abs(diff(v, differences = 2))
  expect_lt(max(d2[(i - 3):(i + 3)]), 10 * median(d2))
  # D1 (validator): just above DBL_MIN, plain-scale pbeta() is off by up to
  # ~10% (pbeta(0.8665, 4975, 25) = 2.2228e-266 vs 2.4013e-266), which gave
  # -611.689 with a 1e-290 threshold; with 1e-240 the Laplace branch is used
  a <- 4975; b <- (1 - 0.995) * 5000
  expect_true(abs(betaregscale:::.brs_obs_loglik(3L, 0.8655, 0.8665, 0.866, a, b) -
                    (-611.615473)) < 1e-4)
  expect_true(abs(cpp1(3L, 0.8655, 0.8665, 0.866, a, b) - (-611.615473)) < 1e-4)
})
