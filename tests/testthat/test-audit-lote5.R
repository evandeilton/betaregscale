# ============================================================================ #
# Lote 5 (2026-09 audit): Armadillo backend, chain-rule derivatives, brsmm
# inner Newton (LM), symmetric-root quadrature, NaN semantics, structural stops
# ============================================================================ #

# brs data for one repar/link combination: delta 0..3 and two far-tail outliers
.l5_brs_data <- function(seed, repar, n = 60L) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n); z1 <- rnorm(n)
  mu <- plogis(0.2 + 0.5 * x1 - 0.3 * x2)
  sh <- brs_repar(mu, if (repar == 2L) 0.15 else 25, repar = if (repar == 0L) 1L else repar)
  Y <- brs_check(round(100 * rbeta(n, sh$shape1, sh$shape2)), ncuts = 100L)
  delta <- as.integer(Y[, "delta"])
  delta[seq(1, n, 10)] <- 0L; delta[seq(2, n, 10)] <- 1L; delta[seq(3, n, 10)] <- 2L
  Y[4:5, "left"] <- c(0.90, 0.02); Y[4:5, "right"] <- Y[4:5, "left"] + 0.01; delta[4:5] <- 3L
  list(X = cbind(1, x1, x2), Z = cbind(1, z1), l = unname(Y[, "left"]),
       r = unname(Y[, "right"]), yt = unname(Y[, "yt"]), delta = delta)
}

# validator's random-slope data (q_re = 2 or 3) as brsmm() builds it
.l5_slope <- function(G, m, seed, q = 2L) {
  set.seed(seed); id <- rep(seq_len(G), each = m); x <- rnorm(G * m); w <- rnorm(G * m)
  b0 <- rnorm(G, 0, .5); b1 <- rnorm(G, 0, .3); b2 <- rnorm(G, 0, .2)
  mu <- plogis(0.2 + (0.5 + b1[id]) * x + (if (q == 3L) (0.3 + b2[id]) * w else 0) + b0[id])
  y <- round(100 * rbeta(G * m, mu * 30, (1 - mu) * 30))
  p <- suppressMessages(suppressWarnings(brs_prep(data.frame(y = y, x = x, w = w, id = factor(id)), ncuts = 100)))
  X <- if (q == 3L) model.matrix(~ x + w, p) else model.matrix(~ x, p)
  list(X = X, Z = matrix(1, nrow(p), 1), Xr = X, l = p$left, r = p$right, yt = p$yt,
       delta = as.integer(p$delta), group = as.integer(factor(p$id)))
}
.l5_ll <- function(P, par, me, np) betaregscale:::.brsmm_loglik_eigen(
  par, P$X, P$Z, P$Xr, P$l, P$r, P$yt, P$delta, P$group, 0L, 0L, 2L, me, np)
.l5_gr <- function(P, par, me, np) betaregscale:::.brsmm_grad_cpp(
  par, P$X, P$Z, P$Xr, P$l, P$r, P$yt, P$delta, P$group, 0L, 0L, 2L, me, np)
.l5_he <- function(P, par, me, np) betaregscale:::.brsmm_hessian_cpp(
  par, P$X, P$Z, P$Xr, P$l, P$r, P$yt, P$delta, P$group, 0L, 0L, 2L, me, np)
# Richardson central differences (steps h max(1,|x|) and half)
.l5_fd <- function(f, x, h0) vapply(seq_along(x), function(j) {
  h <- h0 * max(1, abs(x[j])); e <- replace(numeric(length(x)), j, 1)
  a <- (f(x + h * e) - f(x - h * e)) / (2 * h)
  b <- (f(x + h / 2 * e) - f(x - h / 2 * e)) / h
  (4 * b - a) / 3
}, 1)

test_that("L5-1: brs chain-rule gradient and Hessian match numDeriv (28 repar/link combos)", {
  tab <- list("0" = list(mu = c("log", "sqrt"), phi = c("log", "sqrt")),
              "1" = list(mu = c("logit", "probit", "cauchit", "cloglog"), phi = c("log", "sqrt")),
              "2" = list(mu = c("logit", "probit", "cauchit", "cloglog"),
                         phi = c("logit", "probit", "cauchit", "cloglog")))
  ma <- list(eps = 1e-3, d = 1e-3, zero.tol = Inf, r = 6)   # absolute steps
  k <- 0L
  for (rp in names(tab)) for (lmn in tab[[rp]]$mu) for (lpn in tab[[rp]]$phi) {
    k <- k + 1L; repar <- as.integer(rp); d <- .l5_brs_data(100L + k, repar)
    lm <- betaregscale:::link_to_code(lmn); lp <- betaregscale:::link_to_code(lpn)
    ctr <- c(betaregscale:::apply_link(if (repar == 0L) 3 else 0.55, lmn), 0.3, -0.2,
             betaregscale:::apply_link(if (repar == 2L) 0.15 else if (repar == 1L) 25 else 5, lpn), 0.1)
    for (pt in 1:2) {
      set.seed(1000L * k + pt); p <- ctr + rnorm(5, 0, 0.2)
      fv <- function(b) betaregscale:::.brs_loglik_variable_cpp(b, d$X, d$Z, d$l, d$r, d$yt, d$delta, lm, lp, repar)
      g <- betaregscale:::.brs_grad_variable_cpp(p, d$X, d$Z, d$l, d$r, d$yt, d$delta, lm, lp, repar)
      H <- betaregscale:::.brs_hessian_variable_cpp(p, d$X, d$Z, d$l, d$r, d$yt, d$delta, lm, lp, repar)
      expect_lt(max(abs(g - numDeriv::grad(fv, p, method.args = ma))) / max(abs(g)), 1e-6)
      expect_lt(max(abs(H - numDeriv::hessian(fv, p, method.args = ma))) / max(abs(H)), 1e-5)
      ff <- function(b) betaregscale:::.brs_loglik_fixed_cpp(b, d$X, d$l, d$r, d$yt, d$delta, lm, lp, repar)
      gf <- betaregscale:::.brs_grad_fixed_cpp(p[1:4], d$X, d$l, d$r, d$yt, d$delta, lm, lp, repar)
      Hf <- betaregscale:::.brs_hessian_fixed_cpp(p[1:4], d$X, d$l, d$r, d$yt, d$delta, lm, lp, repar)
      expect_lt(max(abs(gf - numDeriv::grad(ff, p[1:4], method.args = ma))) / max(abs(gf)), 1e-6)
      expect_lt(max(abs(Hf - numDeriv::hessian(ff, p[1:4], method.args = ma))) / max(abs(Hf)), 1e-5)
    }
  }
  expect_equal(k, 28L)
})

test_that("L5-2: brsmm gradient and Hessian (q_re = 1 and 2) match finite differences", {
  skip_on_cran()   # ~10 s of numerical differentiation
  set.seed(1); G <- 30; m <- 8; id <- rep(1:G, each = m); x <- rnorm(G * m)
  mu <- plogis(0.2 + 0.5 * x + rnorm(G, 0, 0.5)[id])
  p1 <- suppressMessages(brs_prep(data.frame(y = round(100 * rbeta(G * m, mu * 30, (1 - mu) * 30)),
                                             x = x, id = factor(id)), ncuts = 100L))
  P1 <- list(X = cbind(1, p1$x), Z = matrix(1, nrow(p1), 1), Xr = matrix(1, nrow(p1), 1),
             l = p1$left, r = p1$right, yt = p1$yt, delta = as.integer(p1$delta),
             group = as.integer(factor(p1$id)))
  P2 <- .l5_slope(30L, 10L, 5L)
  cases <- list(list(P1, c(0.25, 0.5, -3.2, -0.4)), list(P2, c(0.3, 0.52, -3.4, -0.6, 0.03, -1.1)))
  for (cs in cases) for (me in 0:2) {
    P <- cs[[1]]; par <- cs[[2]]; np <- c(11L, 11L, 256L)[me + 1]
    f <- function(b) .l5_ll(P, b, me, np)
    g <- .l5_gr(P, par, me, np)
    expect_lt(max(abs(g - .l5_fd(f, par, 1e-2))) / max(abs(g)), 1e-5)
    if (me < 2L) {
      H <- .l5_he(P, par, me, np)
      Hn <- numDeriv::hessian(f, par, method.args = list(eps = 1e-2, d = 1e-2, zero.tol = Inf, r = 3))
      expect_lt(max(abs(H - Hn)) / max(abs(H)), 1e-3)
      expect_true(all(eigen(-H, symmetric = TRUE, only.values = TRUE)$values > 0))
    }
  }
})

test_that("L5-3: structural errors stop in the compiled functions", {
  d <- .l5_brs_data(1L, 2L)
  fx <- function(p, ...) betaregscale:::.brs_loglik_fixed_cpp(p, d$X, d$l, d$r, d$yt, d$delta, 0L, 0L, 2L, ...)
  expect_true(is.finite(fx(c(0.1, 0.2, 0.1, -1))))
  expect_error(fx(c(0.1, 0.2, -1)), "param must have length 4")
  expect_error(fx(c(0.1, 0.2, 0.1, -1, 5)), "param must have length 4")
  expect_error(betaregscale:::.brs_grad_variable_cpp(c(0, 0, 0, 0, 0), d$X, d$Z, replace(d$l, 2, NA),
               d$r, d$yt, d$delta, 0L, 0L, 2L), "must not contain NA")
  expect_error(betaregscale:::.brs_hessian_fixed_cpp(c(0, 0, 0, 0), replace(d$X, 3, Inf), d$l, d$r,
               d$yt, d$delta, 0L, 0L, 2L), "must not contain NA")
  P <- .l5_slope(10L, 6L, 3L); par <- c(0.2, 0.5, -3, log(0.5), 0, log(0.3))
  expect_true(is.finite(.l5_ll(P, par, 0L, 11L)))
  expect_error(.l5_ll(P, par[-1], 0L, 11L), "param must have length 6")
  expect_error(.l5_gr(P, c(par, 1), 0L, 11L), "param must have length 6")
  Pn <- P; Pn$Xr[1, 2] <- NA
  expect_error(.l5_ll(Pn, par, 0L, 11L), "must not contain NA")
  Pg <- P; Pg$group[1] <- 1e9L
  expect_error(.l5_ll(Pg, par, 0L, 11L), "exceeds the number of rows")
})

test_that("L5-4: brsmm SEs are stable at two practically identical optima (validator parE/parA)", {
  skip_on_cran()   # ~10 s of numerical differentiation
  P <- .l5_slope(30L, 10L, 5L)
  parE <- c(0.30975346835122353, 0.52528622453928187, -3.4418635664742023,
            -0.5917219100895259, 0.021629785845768721, -1.0961160956307257)
  parA <- c(0.3097055261782567, 0.52536935194325396, -3.441873721736203,
            -0.59173301576275028, 0.021651711300330344, -1.0960664302822358)
  for (me in 0:1) {
    seE <- sqrt(diag(solve(-.l5_he(P, parE, me, 11L))))
    seA <- sqrt(diag(solve(-.l5_he(P, parA, me, 11L))))
    expect_true(all(is.finite(seE)) && all(is.finite(seA)))
    expect_lt(max(abs(seE - seA) / seE), 1e-3)
  }
})

test_that("L5-5: QMC (q_re = 3) is smooth: analytic gradient matches small-step differences", {
  P <- .l5_slope(25L, 12L, 6L, q = 3L)
  par <- c(0.2, 0.5, 0.3, -3.4, log(0.5), 0.2, -0.1, log(0.3), 0.1, log(0.2)) + 0.1
  f <- function(b) .l5_ll(P, b, 2L, 256L); g <- .l5_gr(P, par, 2L, 256L)
  for (h in c(1e-3, 1e-4)) {
    fd <- vapply(seq_along(par), function(j) { e <- replace(numeric(10), j, h); (f(par + e) - f(par - e)) / (2 * h) }, 1)
    expect_lt(max(abs(fd - g)) / max(abs(g)), 1e-4)
  }
})

test_that("L5-6: LM inner modes are stationary with positive-definite curvature (q_re = 3)", {
  P <- .l5_slope(25L, 12L, 6L, q = 3L)
  base <- c(0.2, 0.5, 0.3, -3.4, log(0.5), 0.2, -0.1, log(0.3), 0.1, log(0.2))
  set.seed(99); pts <- c(list(base), lapply(1:3, function(k) base + rnorm(10, 0, 0.3)))
  for (par in pts) for (warm in c(FALSE, TRUE)) {
    dg <- betaregscale:::.brsmm_mode_diag_cpp(par, P$X, P$Z, P$Xr, P$l, P$r, P$yt, P$delta,
                                              P$group, 0L, 0L, 2L, warm)
    expect_true(all(dg$ok == 1))
    expect_lt(max(dg$grad_inf), 1e-8)
    expect_true(all(dg$min_eig > 0))
  }
})

test_that("L5-7: NaN parameters give the penalty; +-Inf map to the bounds (C++ = R mirror)", {
  set.seed(1); n <- 20L; l <- runif(n, 0.1, 0.8); r <- l + 0.01; yt <- l + 0.005; d <- rep(3L, n)
  X <- matrix(1, n, 1); code <- betaregscale:::link_to_code
  f <- function(p, lm, lp, rp) betaregscale:::.brs_loglik_fixed_cpp(p, X, l, r, yt, d, lm, lp, rp)
  expect_equal(f(c(NaN, 0), code("logit"), code("logit"), 2L), -1e6 * n)
  expect_equal(f(c(0, NaN), code("logit"), code("log"), 1L), -1e6 * n)
  expect_equal(f(c(NaN, 1), code("sqrt"), code("sqrt"), 0L), -1e6 * n)   # sqrt link propagates NaN
  # identity links: param is (mu, phi) itself, so C++ clamps can be compared with the R mirror
  for (rp in 0:2) for (mu in c(-Inf, Inf, 0.3)) for (phi in c(-Inf, Inf, 0.2)) {
    m <- betaregscale:::.clamp_mu_by_repar(mu, rp); ph <- betaregscale:::.clamp_phi_by_repar(phi, rp)
    sh <- switch(as.character(rp), "0" = c(m, ph), "1" = c(m * ph, (1 - m) * ph),
                 "2" = c(m * (1 - ph) / ph, (1 - m) * (1 - ph) / ph))
    sh <- pmin(pmax(sh, 1e-12), 1e8)
    ref <- sum(betaregscale:::.brs_obs_loglik(d, l, r, yt, sh[1], sh[2]))
    expect_equal(f(c(mu, phi), code("identity"), code("identity"), rp), ref, tolerance = 1e-10)
  }
  expect_true(is.nan(betaregscale:::.clamp_mu_by_repar(NaN, 2L)))
  expect_true(is.na(betaregscale:::.clamp_phi_by_repar(NA_real_, 1L)))
  expect_equal(betaregscale:::.clamp_mu_by_repar(-Inf, 0L), 1e-5)
})

test_that("L5-8: AGHQ converges to brute-force integration as n_points grows (q_re = 2)", {
  P <- .l5_slope(4L, 10L, 5L); par <- c(0.2, 0.5, -3.4, log(0.5), 0.2, log(0.3))
  modes <- betaregscale:::.brsmm_group_modes_eigen(par, P$X, P$Z, P$Xr, P$l, P$r, P$yt, P$delta,
                                                   P$group, 0L, 0L, 2L)
  L <- matrix(c(exp(par[4]), par[5], 0, exp(par[6])), 2); D <- L %*% t(L); Pm <- solve(D)
  u <- seq(-9, 9, by = 0.15); U <- as.matrix(expand.grid(u, u)); bf <- 0
  for (g in 1:4) {
    i <- which(P$group == g); K <- nrow(U); B <- sweep(U * 0.5, 2, modes[g, ], "+")
    eta <- rep(drop(P$X[i, ] %*% par[1:2]), each = K) + as.vector(B %*% t(P$Xr[i, ]))
    mu <- pmin(pmax(plogis(eta), 1e-5), 1 - 1e-5); ph <- plogis(par[3])
    lv <- betaregscale:::.brs_obs_loglik(rep(P$delta[i], each = K), rep(P$l[i], each = K),
            rep(P$r[i], each = K), rep(P$yt[i], each = K), mu * (1 - ph) / ph, (1 - mu) * (1 - ph) / ph)
    h <- rowSums(matrix(lv, K)) - 0.5 * (2 * log(2 * pi) + log(det(D)) + rowSums((B %*% Pm) * B))
    bf <- bf + max(h) + log(sum(exp(h - max(h)))) + 2 * log(0.15 * 0.5)
  }
  e <- vapply(c(3L, 7L, 15L), function(np) abs(.l5_ll(P, par, 1L, np) - bf), 1)
  expect_true(e[2] < e[1] && e[3] < e[2])
  expect_lt(e[3], 1e-6)
})

test_that("L5-9: a brsmm fit does not depend on earlier fits in the session (warm-start reset)", {
  mk <- function(seed) {
    set.seed(seed); G <- 12; m <- 8; id <- rep(seq_len(G), each = m); x <- rnorm(G * m)
    mu <- plogis(0.2 + 0.5 * x + rnorm(G, 0, 0.5)[id])
    y <- round(100 * rbeta(G * m, mu * 30, (1 - mu) * 30))
    suppressMessages(brs_prep(data.frame(y = y, x = x, id = factor(id)), ncuts = 100L))
  }
  dA <- mk(1L); dB <- mk(2L)   # same G and q_re, different data
  fit <- function(d) suppressWarnings(brsmm(y ~ x, random = ~ 1 | id, data = d))
  a1 <- fit(dA); invisible(fit(dB)); a2 <- fit(dA)
  expect_identical(a1$par, a2$par)
  expect_identical(as.numeric(logLik(a1)), as.numeric(logLik(a2)))
  expect_identical(vcov(a1), vcov(a2))
})

test_that("L5-10: NA in delta stops in R and in the compiled code (no NaN -> integer UB)", {
  set.seed(2); d <- data.frame(y = round(100 * rbeta(40, 2, 3)), x = rnorm(40))
  p <- suppressMessages(brs_prep(d, ncuts = 100L))
  p$delta[3] <- NA
  expect_error(brs(y ~ x, data = p), "delta.*NA or other value.*first: 3")
  p$delta[3] <- 7L
  expect_error(brs(y ~ x, data = p), "delta.*must be 0, 1, 2 or 3")
  q <- suppressMessages(brs_prep(d, ncuts = 100L)); X <- cbind(1, q$x)
  dl <- as.integer(q$delta); dl[3] <- NA
  for (f in list(betaregscale:::.brs_loglik_fixed_cpp, betaregscale:::.brs_grad_fixed_cpp,
                 betaregscale:::.brs_hessian_fixed_cpp)) {
    expect_error(f(c(0, 0, 0), X, q$left, q$right, q$yt, dl, 0L, 0L, 2L), "found NA at row 3")
  }
  Z <- matrix(1, nrow(q), 1)
  for (f in list(betaregscale:::.brs_loglik_variable_cpp, betaregscale:::.brs_grad_variable_cpp,
                 betaregscale:::.brs_hessian_variable_cpp)) {
    expect_error(f(c(0, 0, 0), X, Z, q$left, q$right, q$yt, dl, 0L, 0L, 2L), "found NA at row 3")
  }
  # a double delta with NaN is coerced by R (NaN -> NA), never cast in C++
  expect_error(betaregscale:::.brs_loglik_fixed_cpp(c(0, 0, 0), X, q$left, q$right, q$yt,
               replace(as.numeric(q$delta), 3, NaN), 0L, 0L, 2L), "found NA at row 3")
  g <- rep(1:4, each = 10)
  expect_error(betaregscale:::.brsmm_loglik_eigen(c(0, 0, 0, 0), X, Z, Z, q$left, q$right, q$yt,
               dl, g, 0L, 0L, 2L, 0L, 11L), "delta.*found NA at row 3")
})
