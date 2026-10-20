# Clustered simulation, simulate() methods, mixed/cell residuals, BIC(n), folds

# Design of G groups with nj visits and centered time
sim_design <- function(G, nj) {
  d <- expand.grid(visit = seq_len(nj), id = seq_len(G))
  d$tc <- (d$visit - (nj + 1) / 2) / 2
  d
}

# Random-intercept data and its fit: G = 200 groups of 5 rows
sim_fit_ri <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      set.seed(3)
      s <- brsmm_sim(~ tc, random = ~ 1 | id, data = sim_design(200, 5),
                     beta = c(-0.5, 0.4), phi = qlogis(0.2), re_sd = 0.8, ncuts = 10)
      cache <<- list(s = s, m = brsmm(y ~ tc, random = ~ 1 | id, data = s))
    }
    cache
  }
})

test_that("brsmm_sim is reproducible and follows the documented draw order", {
  d <- sim_design(30, 4)
  args <- list(~ tc, random = ~ 1 | id, data = d, beta = c(-0.3, 0.5),
               phi = qlogis(0.25), re_sd = 0.7, ncuts = 10, cutpoints = c(0.08, 0.93))
  set.seed(11)
  s1 <- do.call(brsmm_sim, args)
  set.seed(11)
  s2 <- do.call(brsmm_sim, args)
  expect_identical(s1, s2)
  # By hand: G normals for b, then n beta values, coarsened with the cutpoints
  set.seed(11)
  b <- 0.7 * rnorm(30)
  mu <- plogis(-0.3 + 0.5 * d$tc + b[d$id])
  sh <- brs_repar(mu = mu, phi = 0.25, repar = 2)
  y <- rbeta(nrow(d), sh$shape1, sh$shape2)
  s <- pmin(pmax(round(10 * y), 1), 9)
  s[y < 0.08] <- 0
  s[y >= 0.93] <- 10
  expect_equal(s1$y, s)
  expect_equal(as.numeric(attr(s1, "ranef")), b)
  expect_equal(rownames(attr(s1, "ranef")), as.character(1:30))
  expect_identical(attr(s1, "ncuts"), 10L)
  expect_true(isTRUE(attr(s1, "is_prepared")))
  expect_equal(names(s1)[1:5], c("left", "right", "yt", "y", "delta"))
  expect_true(all(c("visit", "id", "tc") %in% names(s1)))
  # Cells and censoring of the map, as brs_sim() returns for score data
  chk <- .brs_check_core(s1$y, 10, 0.5, NULL, "mid")
  expect_equal(as.matrix(s1[, 1:5]), chk, ignore_attr = TRUE)
  # True cutpoint parameters map back to the cutpoints
  tr <- attr(s1, "truth")
  g <- tr$par[c("(cut_lo)", "(cut_hi)")]
  expect_equal(unname(.brs_cut_values(g, 10, .brs_cut_anchors(10, "mid"))$value),
               c(0.08, 0.93))
})

test_that("simulated random effects have the requested covariance", {
  G <- 20000
  d <- data.frame(id = seq_len(G), tc = rnorm(G))
  sd0 <- c(0.8, 0.3)
  rho <- 0.3
  set.seed(5)
  s <- brsmm_sim(~ tc, random = ~ 1 + tc | id, data = d, beta = c(0, 0.2),
                 phi = qlogis(0.2), re_sd = sd0, re_cor = rho, ncuts = 10)
  b <- attr(s, "ranef")
  D <- attr(s, "truth")$D
  expect_equal(unname(D), matrix(c(0.64, 0.072, 0.072, 0.09), 2), tolerance = 1e-12)
  V <- stats::cov(b)
  # Monte Carlo SE of a sample (co)variance: sqrt((D_ii D_jj + D_ij^2) / G)
  se <- sqrt((outer(diag(D), diag(D)) + D^2) / G)
  expect_true(all(abs(V - D) < 4 * se))
  expect_equal(colnames(b), c("(Intercept)", "tc"))
})

test_that("simulated score frequencies match the marginal pmf of brs_marginal()", {
  # Fit only to get a brsmm object; its parameters are set to the true ones
  set.seed(1)
  small <- brsmm_sim(~ tc, random = ~ 1 | id, data = sim_design(40, 5),
                     beta = c(-0.5, 0.4), phi = qlogis(0.2), re_sd = 0.9, ncuts = 6)
  m <- brsmm(y ~ tc, random = ~ 1 | id, data = small)
  tr <- attr(small, "truth")$par
  m$par[names(tr)] <- tr
  pm <- brs_marginal(m, newdata = data.frame(tc = 0.5), type = "prob")$estimate
  expect_equal(sum(pm), 1, tolerance = 1e-10)
  # One row per group: independent draws from the marginal law
  G <- 40000
  set.seed(2)
  big <- brsmm_sim(~ tc, random = ~ 1 | id, data = data.frame(id = seq_len(G), tc = 0.5),
                   beta = c(-0.5, 0.4), phi = qlogis(0.2), re_sd = 0.9, ncuts = 6)
  obs <- tabulate(big$y + 1, 7)
  expect_gt(stats::chisq.test(obs, p = pm)$p.value, 1e-3)
  # Mean score within 4 Monte Carlo SEs
  es <- sum(0:6 * pm)
  sds <- sqrt(sum((0:6 - es)^2 * pm))
  expect_lt(abs(mean(big$y) - es), 4 * sds / sqrt(G))
})

test_that("invalid random-effect specifications are errors", {
  d <- sim_design(10, 3)
  d$x2 <- rnorm(nrow(d))
  base <- list(~ tc, data = d, beta = c(0, 0.1), phi = qlogis(0.2), ncuts = 10)
  call_sim <- function(...) do.call(brsmm_sim, c(base, list(...)))
  expect_error(call_sim(random = ~ 1 | id, re_sd = -1), "re_sd")
  expect_error(call_sim(random = ~ 1 | id, re_sd = 0), "re_sd")
  expect_error(call_sim(random = ~ 1 | id, re_sd = c(1, 1)), "re_sd")
  expect_error(call_sim(random = ~ 1 | id, re_sd = NA), "re_sd")
  expect_error(call_sim(random = ~ 1 | id), "re_sd")
  expect_error(call_sim(random = ~ 1 + tc | id, re_sd = c(1, 0.5), re_cor = 1.2),
               "re_cor")
  expect_error(call_sim(random = ~ 1 + tc | id, re_sd = c(1, 0.5),
                        re_cor = matrix(c(1, 0.2, 0.3, 1), 2)), "re_cor")
  expect_error(call_sim(random = ~ 1 + tc | id, re_sd = c(1, 0.5), re_cor = 1),
               "positive definite")
  R <- matrix(c(1, 0.9, -0.9, 0.9, 1, 0.9, -0.9, 0.9, 1), 3)
  expect_error(do.call(brsmm_sim, c(base, list(random = ~ 1 + tc + x2 | id,
                                               re_sd = c(1, 1, 1), re_cor = R))),
               "positive definite")
  expect_error(call_sim(random = ~ 1 | grp, re_sd = 1), "not found")
  expect_error(call_sim(random = ~ 1 | id, re_sd = 1, cutpoints = c(0.5, 0.9)),
               "cutpoints")
})

test_that("brsmm() recovers the parameters of brsmm_sim() (G = 200, n_j = 5)", {
  fx <- sim_fit_ri()
  tr <- attr(fx$s, "truth")$par
  expect_identical(names(tr), names(coef(fx$m)))
  z <- (coef(fx$m) - tr) / sqrt(diag(vcov(fx$m)))
  expect_true(all(abs(z) < 4))
})

test_that("simulate() for brsmm matches brsmm_sim() and conditions on the modes", {
  fx <- sim_fit_ri()
  m <- fx$m
  p <- m$par
  set.seed(9)
  s_new <- brsmm_sim(~ tc, random = ~ 1 | id, data = fx$s,
                     beta = unname(p[1:2]), phi = unname(p[3]),
                     re_sd = unname(m$random$sd_b), ncuts = 10)
  sim <- simulate(m, nsim = 2, seed = 9)
  expect_equal(sim$sim_1, s_new$y)
  expect_equal(attr(sim, "ranef")[[1]], attr(s_new, "ranef"), tolerance = 1e-12)
  expect_identical(attr(sim, "seed")[1], 9)
  # The caller's RNG stream is untouched with a seed
  set.seed(4); a <- runif(1)
  set.seed(4); invisible(simulate(m, seed = 1)); b <- runif(1)
  expect_identical(a, b)
  # Conditional: random effects at the modes, only the beta draws vary
  sc <- simulate(m, nsim = 1, seed = 2, re = "conditional")
  set.seed(2)
  eta <- p[1] + p[2] * fx$s$tc + m$random$mode_b[m$group_index]
  sh <- .brs_shapes_cpp(eta, rep(p[3], length(eta)), m$link, m$link_phi, m$repar)
  expect_equal(sc$sim_1, round(10 * rbeta(length(eta), sh$a, sh$b)))
  expect_null(attr(sc, "ranef"))
  expect_error(simulate(m, nsim = 0), "nsim")
})

test_that("simulate() for brs matches brs_sim() and uses estimated cutpoints", {
  set.seed(31)
  d <- data.frame(x = rnorm(400))
  s <- brs_sim(~ x, data = d, beta = c(0.2, 0.5), phi = qlogis(0.2), ncuts = 10)
  f <- brs(y ~ x, data = s)
  set.seed(6)
  s2 <- brs_sim(~ x, data = s, beta = unname(coef(f)[1:2]), phi = unname(coef(f)[3]),
                ncuts = 10)
  sim <- simulate(f, nsim = 3, seed = 6)
  expect_equal(sim$sim_1, s2$y)
  expect_identical(dim(sim), c(400L, 3L))
  expect_identical(simulate(f, nsim = 3, seed = 6), sim)
  # Boundary fit: the estimated cutpoints coarsen the draws
  set.seed(32)
  sb <- brs_sim(~ x, data = d, beta = c(0.2, 0.5), phi = qlogis(0.2), ncuts = 10,
                cutpoints = c(0.1, 0.88))
  fb <- suppressWarnings(brs(y ~ x, data = sb, cutpoints = "boundary"))
  set.seed(7)
  s3 <- brs_sim(~ x, data = sb, beta = unname(coef(fb)[1:2]), phi = unname(coef(fb)[3]),
                ncuts = 10, cutpoints = fb$cut$value)
  expect_equal(simulate(fb, seed = 7)$sim_1, s3$y)
  # Exact rows keep the latent value
  se <- brs_sim(~ x, data = d[1:100, , drop = FALSE], beta = c(0.2, 0.5),
                phi = qlogis(0.2), ncuts = 10, delta = 0)
  fe <- brs(y ~ x, data = se)
  ye <- simulate(fe, seed = 1)$sim_1
  expect_true(all(ye > 0 & ye < 1 & ye != round(ye)))
})

test_that("marginal and posterior quantile residuals are close to N(0, 1)", {
  fx <- sim_fit_ri()
  m <- fx$m
  set.seed(20)
  rm <- residuals(m, type = "rqr", re = "marginal")
  rp <- residuals(m, type = "rqr", re = "posterior")
  for (r in list(rm, rp)) {
    expect_length(r, nrow(fx$s))
    expect_true(all(is.finite(r)))
    expect_lt(abs(mean(r)), 0.1)
    expect_lt(abs(stats::sd(r) - 1), 0.08)
    expect_gt(stats::ks.test(r, "pnorm")$p.value, 1e-3)
  }
  # Reproducible with a seed; the default is the conditional residual
  set.seed(21); a <- residuals(m, type = "rqr", re = "posterior")
  set.seed(21); b <- residuals(m, type = "rqr", re = "posterior")
  expect_identical(a, b)
  set.seed(22); a <- residuals(m, type = "rqr")
  set.seed(22); b <- residuals(m, type = "rqr", re = "conditional")
  expect_identical(a, b)
  expect_error(residuals(m, type = "pearson", re = "marginal"), "applies to")
  expect_error(residuals(m, type = "rqr", re = "posterior", n_draws = 0), "n_draws")
})

test_that("marginal quantile residuals use the marginal pmf of the score", {
  fx <- sim_fit_ri()
  m <- fx$m
  rows <- c(1, 7, 13)
  nd <- fx$s[rows, ]
  pm <- matrix(brs_marginal(m, newdata = nd, type = "prob")$estimate, ncol = 11,
               byrow = TRUE)
  cdf <- t(apply(pm, 1L, cumsum))
  s <- nd$y
  lo <- ifelse(s == 0, 0, cdf[cbind(seq_along(s), pmax(s, 1))])
  hi <- cdf[cbind(seq_along(s), s + 1)]
  set.seed(23)
  r <- residuals(m, type = "rqr", re = "marginal")[rows]
  u <- stats::pnorm(r)
  expect_true(all(u >= lo - 1e-6 & u <= hi + 1e-6))
  # Cell deviance from the same marginal probabilities
  dm <- residuals(m, type = "cell_deviance", re = "marginal")[rows]
  expect_equal(unname(dm^2), -2 * log(pm[cbind(seq_along(s), s + 1)]),
               tolerance = 1e-6)
})

test_that("posterior draws of a random intercept follow its posterior", {
  set.seed(41)
  s <- brsmm_sim(~ 1, random = ~ 1 | id, data = sim_design(15, 3), beta = -0.3,
                 phi = qlogis(0.3), re_sd = 0.8, ncuts = 6)
  m <- brsmm(y ~ 1, random = ~ 1 | id, data = s)
  et <- .brs_res_eta(m)
  # Posterior of group 1 on a fine grid
  i1 <- which(m$group_index == 1L)
  bg <- seq(-6, 6, length.out = 4001)
  lpost <- vapply(bg, function(b) {
    sh <- .brs_shapes_cpp(et$mu[i1] + b, et$phi[i1], m$link, m$link_phi, m$repar)
    sum(.brs_obs_loglik(m$Y[i1, "delta"], m$Y[i1, "left"], m$Y[i1, "right"],
                        m$Y[i1, "yt"], sh$a, sh$b))
  }, numeric(1)) + stats::dnorm(bg, 0, m$random$sd_b, log = TRUE)
  w <- exp(lpost - max(lpost))
  w <- w / sum(w)
  pmean <- sum(w * bg)
  psd <- sqrt(sum(w * (bg - pmean)^2))
  set.seed(42)
  draws <- replicate(400, .brsmm_posterior_b(m, et$mu, et$phi, 200L)[1, 1])
  expect_lt(abs(mean(draws) - pmean), 4 * psd / sqrt(400))
  expect_lt(abs(stats::sd(draws) / psd - 1), 0.15)
})

test_that("cell deviance residuals add up to -2 logLik and carry the score sign", {
  set.seed(51)
  d <- data.frame(x = rnorm(300), z = runif(300))
  s <- brs_sim(~ x | z, data = d, beta = c(0.2, 0.6), zeta = c(-1.5, 0.5), ncuts = 10)
  f <- brs(y ~ x | z, data = s)
  r <- residuals(f, type = "cell_deviance")
  expect_equal(sum(r^2), -2 * as.numeric(logLik(f)), tolerance = 1e-8)
  es <- brs_marginal(f, type = "expected_score")$estimate
  expect_identical(unname(sign(r)), ifelse(s$y >= es, 1, -1))
  # Right map and boundary cutpoints
  sr <- brs_sim(~ x, data = d, beta = c(0.2, 0.6), phi = qlogis(0.2), ncuts = 10,
                interval = "right")
  fr <- brs(y ~ x, data = sr)
  expect_equal(sum(residuals(fr, type = "cell_deviance")^2),
               -2 * as.numeric(logLik(fr)), tolerance = 1e-8)
  sb <- brs_sim(~ x, data = d, beta = c(0.2, 0.6), phi = qlogis(0.2), ncuts = 10,
                cutpoints = c(0.12, 0.9))
  fb <- suppressWarnings(brs(y ~ x, data = sb, cutpoints = "boundary"))
  expect_equal(sum(residuals(fb, type = "cell_deviance")^2),
               -2 * as.numeric(logLik(fb)), tolerance = 1e-8)
  # Exact observations have no cell: NA
  se <- brs_sim(~ x, data = d[1:80, , drop = FALSE], beta = c(0.2, 0.6),
                phi = qlogis(0.2), ncuts = 10, delta = 0)
  expect_true(all(is.na(residuals(brs(y ~ x, data = se), type = "cell_deviance"))))
  # brsmm, conditional: -2 log P(s | modes)
  fx <- sim_fit_ri()
  m <- fx$m
  dc <- residuals(m, type = "cell_deviance")
  et <- .brs_res_eta(m)
  sh <- .brs_shapes_cpp(et$mu + m$random$mode_b[m$group_index], et$phi, m$link,
                        m$link_phi, m$repar)
  lp <- .brs_obs_loglik(m$Y[, "delta"], m$Y[, "left"], m$Y[, "right"], m$Y[, "yt"],
                        sh$a, sh$b)
  expect_equal(unname(dc^2), -2 * lp)
})

test_that("BIC(n = 'groups') with a brs fit stops in either order", {
  fx <- sim_fit_ri()
  m <- fx$m
  f <- brs(y ~ tc, data = fx$s)
  expect_error(BIC(f, m, n = "groups"), "needs \"brsmm\" fits only")
  expect_error(BIC(m, f, n = "groups"), "needs \"brsmm\" fits only")
  expect_error(BIC(f, n = "groups"), "needs \"brsmm\" fits only")
  expect_identical(BIC(f, m, n = "obs"), BIC(f, m))
  expect_identical(rownames(BIC(f, m, n = "obs")), c("f", "m"))
  expect_identical(BIC(f, n = "obs"), BIC(f))
})

test_that("user folds keep their type whether given by row or by group", {
  set.seed(61)
  s <- brsmm_sim(~ 1, random = ~ 1 | id, data = sim_design(12, 2), beta = -0.3,
                 phi = qlogis(0.3), re_sd = 0.5, ncuts = 6)
  grp <- factor(s$id)
  by_row_num <- .brsmm_cv_user_folds(rep(1:3, length.out = 12)[s$id], grp, FALSE)
  by_name_num <- .brsmm_cv_user_folds(stats::setNames(rep(1:3, length.out = 12),
                                                      1:12), grp, TRUE)
  expect_identical(by_row_num, by_name_num)
  expect_type(by_row_num, "integer")
  lab <- c("a", "b", "c")[rep(1:3, length.out = 12)]
  by_row_chr <- .brsmm_cv_user_folds(lab[s$id], grp, FALSE)
  expect_identical(by_row_chr, .brsmm_cv_user_folds(stats::setNames(lab, 1:12), grp, TRUE))
  expect_identical(.brsmm_cv_user_folds(factor(lab[s$id]), grp, FALSE), by_row_chr)
  expect_type(by_row_chr, "character")
  expect_identical(names(by_row_chr), levels(grp))
})

test_that("posterior curvature at the modes matches numDeriv", {
  set.seed(71)
  s <- brsmm_sim(~ tc, random = ~ 1 + tc | id, data = sim_design(30, 5),
                 beta = c(-0.4, 0.3), phi = qlogis(0.2), re_sd = c(0.8, 0.5),
                 ncuts = 10)
  m <- suppressWarnings(brsmm(y ~ tc, random = ~ 1 + tc | id, data = s))
  et <- .brs_res_eta(m)
  lp <- .brsmm_logpost_fun(m, et$mu, et$phi)
  B0 <- m$random$mode_b
  H <- .brsmm_curvature(lp, B0)
  for (g in c(1L, 17L)) {
    f <- function(b) {
      B <- B0
      B[g, ] <- b
      lp(B)[g]
    }
    expect_equal(-H[g, , ], numDeriv::hessian(f, B0[g, ]), tolerance = 1e-5)
    # The modes are stationary points of the log posterior
    expect_lt(max(abs(numDeriv::grad(f, B0[g, ]))), 1e-4)
  }
})

test_that("folds given as a factor named by group level are read by name", {
  set.seed(62)
  s <- brsmm_sim(~ 1, random = ~ 1 | id, data = sim_design(12, 2), beta = -0.3,
                 phi = qlogis(0.3), re_sd = 0.5, ncuts = 6)
  grp <- factor(s$id)
  lab <- c("a", "b", "c")[rep(1:3, length.out = 12)]
  ff <- stats::setNames(factor(lab), rev(levels(grp)))
  got <- .brsmm_cv_user_folds(ff, grp, by_name = TRUE)
  expect_identical(got, stats::setNames(lab, rev(levels(grp)))[levels(grp)])
})

test_that("posterior residuals work with one proposal; marginal pieces by chunks", {
  fx <- sim_fit_ri()
  m <- fx$m
  set.seed(24)
  r1 <- residuals(m, type = "rqr", re = "posterior", n_draws = 1)
  expect_length(r1, nrow(fx$s))
  expect_true(all(is.finite(r1)))
  et <- .brs_res_eta(m)
  one <- .brsmm_res_marginal_parts(m, et$mu, et$phi, need_es = TRUE)
  many <- .brsmm_res_marginal_parts(m, et$mu, et$phi, need_es = TRUE,
                                    max_cells = 5000)
  expect_equal(one, many, tolerance = 1e-12)
})
