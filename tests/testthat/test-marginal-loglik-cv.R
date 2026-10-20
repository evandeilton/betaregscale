# Marginal predictions, log-likelihood contributions, group CV and BIC(n)

# NRS-11-like scores: `G` patients x 3 times, random intercept (and slope on tc)
mm_test_data <- function(G = 40, sd0 = 0.8, sd1 = 0, seed = 21) {
  set.seed(seed)
  d <- expand.grid(id = seq_len(G), time = c("6h", "12h", "24h"))
  d$tc <- c(-1, 0, 1)[d$time]
  b0 <- rnorm(G, sd = sd0)
  b1 <- rnorm(G, sd = sd1)
  eta <- -1.3 + c(0, 0.75, 0.3)[d$time] + b0[d$id] + b1[d$id] * d$tc
  shp <- brs_repar(mu = plogis(eta), phi = 0.2, repar = 2)
  d$y <- round(10 * rbeta(nrow(d), shp$shape1, shp$shape2))
  d
}

nd3 <- data.frame(time = c("6h", "12h", "24h"), tc = c(-1, 0, 1))

# Conditional cell probabilities of one row at linear predictor values `eta`
cond_pmf <- function(fit, eta, eta_phi, theta = fit$par) {
  e <- .brs_marg_ends(fit, theta)
  sh <- .brs_shapes_cpp(eta, rep_len(eta_phi, length(eta)), fit$link,
                        fit$link_phi, fit$repar)
  .brs_cell_probs(sh$a, sh$b, e$left, e$right, fit$ncuts)
}

test_that("marginal pmf rows sum to 1 and agree with brute-force integration (q = 1)", {
  d <- mm_test_data()
  m <- brsmm(y ~ time, random = ~ 1 | id, data = d, ncuts = 10,
             int_method = "aghq", n_points = 15)
  pm <- brs_marginal(m, nd3)
  expect_s3_class(pm, "brs_marginal")
  P <- matrix(pm$estimate, nrow = 3, byrow = TRUE)
  expect_equal(rowSums(P), rep(1, 3), tolerance = 1e-12)
  expect_true(all(P >= 0))

  beta <- coef(m, "mean")
  sdb <- exp(coef(m, "random"))
  eta_phi <- coef(m, "precision")
  X <- model.matrix(~ time, data.frame(time = factor(nd3$time, levels = c("6h", "12h", "24h"))))
  for (i in 1:3) {
    eta0 <- sum(X[i, ] * beta)
    ex <- vapply(0:10, function(s) {
      stats::integrate(function(b) {
        cond_pmf(m, eta0 + b, eta_phi)[, s + 1L] * dnorm(b, 0, sdb)
      }, -Inf, Inf, rel.tol = 1e-12, abs.tol = 1e-14)$value
    }, numeric(1))
    expect_equal(P[i, ], ex, tolerance = 1e-9)
  }
})

test_that("marginal pmf with two random effects matches 2-D quadrature and Monte Carlo", {
  d <- mm_test_data(G = 50, sd1 = 0.5, seed = 4)
  m <- suppressWarnings(brsmm(y ~ time, random = ~ 1 + tc | id, data = d, ncuts = 10))
  pm <- brs_marginal(m, nd3)
  P <- matrix(pm$estimate, nrow = 3, byrow = TRUE)
  expect_equal(rowSums(P), rep(1, 3), tolerance = 1e-12)

  L <- m$random$L
  beta <- coef(m, "mean")
  eta_phi <- coef(m, "precision")
  X <- model.matrix(~ time, data.frame(time = factor(nd3$time, levels = c("6h", "12h", "24h"))))
  Xr <- cbind(1, nd3$tc)
  # Product Gauss-Hermite rule in b = L z (independent of the 1-D reduction)
  gh <- .brs_gh_normal(40)
  zz <- as.matrix(expand.grid(gh$x, gh$x))
  ww <- as.vector(outer(gh$w, gh$w))
  bb <- zz %*% t(L)
  set.seed(99)
  bmc <- matrix(rnorm(2e5 * 2), ncol = 2) %*% t(L)
  for (i in 1:3) {
    eta0 <- sum(X[i, ] * beta)
    P2 <- colSums(ww * cond_pmf(m, eta0 + drop(bb %*% Xr[i, ]), eta_phi))
    expect_equal(P[i, ], P2, tolerance = 1e-8)
    Pmc <- cond_pmf(m, eta0 + drop(bmc %*% Xr[i, ]), eta_phi)
    mc <- colMeans(Pmc)
    mc_se <- apply(Pmc, 2, sd) / sqrt(nrow(Pmc))
    expect_true(all(abs(P[i, ] - mc) < 5 * mc_se + 1e-12))
  }
})

test_that("tiny random-effect SD gives the conditional values; brs marginal = conditional", {
  d <- mm_test_data()
  m <- brsmm(y ~ time, random = ~ 1 | id, data = d, ncuts = 10)
  m0 <- m
  m0$par[grep("re_chol", names(m0$par))] <- -12
  for (tp in c("prob", "expected_score", "exceed")) {
    a <- brs_marginal(m0, nd3, type = tp)$estimate
    b <- brs_marginal(m0, nd3, type = tp, conditional = TRUE)$estimate
    expect_equal(a, b, tolerance = 1e-9)
  }
  # With the fitted SD they differ (nonlinear link)
  a <- brs_marginal(m, nd3, type = "expected_score")$estimate
  b <- brs_marginal(m, nd3, type = "expected_score", conditional = TRUE)$estimate
  expect_true(all(abs(a - b) > 1e-3))
  # Conditional values at b = 0 equal predict() for a new patient
  expect_equal(b, predict(m, newdata = transform(nd3, id = 999), type = "expected_score"),
               tolerance = 1e-10)

  f <- brs(y ~ time, data = d, ncuts = 10)
  e1 <- brs_marginal(f, nd3, type = "expected_score")
  e2 <- brs_marginal(f, nd3, type = "expected_score", conditional = TRUE)
  expect_identical(e1$estimate, e2$estimate)
  expect_equal(e1$estimate, predict(f, newdata = nd3, type = "expected_score"),
               tolerance = 1e-9)
  pf <- brs_marginal(f, nd3)
  expect_equal(matrix(pf$estimate, 3, byrow = TRUE),
               unname(brs_predict_scoreprob(f, newdata = nd3)), tolerance = 1e-9,
               ignore_attr = TRUE)
})

test_that("expected score, exceedance and contrasts follow from the pmf", {
  d <- mm_test_data()
  m <- brsmm(y ~ time, random = ~ 1 | id, data = d, ncuts = 10)
  P <- matrix(brs_marginal(m, nd3)$estimate, nrow = 3, byrow = TRUE)
  es <- brs_marginal(m, nd3, type = "expected_score")
  expect_equal(es$estimate, drop(P %*% 0:10), tolerance = 1e-12)
  ex <- brs_marginal(m, nd3, type = "exceed", k = c(0, 4, 10))
  expect_equal(ex$estimate, as.vector(t(cbind(rowSums(P), rowSums(P[, 5:11]), P[, 11]))),
               tolerance = 1e-12)
  expect_equal(ex$std.error[ex$k == 0], rep(0, 3), tolerance = 1e-8)

  cs <- rbind(c(2, 1), c(3, 1))
  ct <- brs_marginal(m, nd3, type = "expected_score", contrast = cs)
  expect_equal(ct$estimate, es$estimate[c(2, 3)] - es$estimate[1], tolerance = 1e-12)
  J <- attr(es, "jacobian")
  V <- vcov(m)
  Jd <- J[c(2, 3), ] - J[c(1, 1), ]
  expect_equal(ct$std.error, sqrt(diag(Jd %*% V %*% t(Jd))), tolerance = 1e-8)
  expect_equal(ct$lower, ct$estimate - qnorm(0.975) * ct$std.error)
  # The three contrast forms agree
  c2 <- brs_marginal(m, nd3, type = "expected_score", contrast = list(c(2, 1), c(3, 1)))
  expect_equal(c2$estimate, ct$estimate)
  c3 <- brs_marginal(m, nd3, type = "expected_score", contrast = c(2, 1))
  expect_equal(c3$estimate, ct$estimate[1])
  ctp <- brs_marginal(m, nd3, type = "exceed", k = 4, contrast = cs)
  expect_equal(ctp$estimate, ex$estimate[ex$k == 4][2:3] - ex$estimate[ex$k == 4][1],
               tolerance = 1e-12)
  expect_error(brs_marginal(m, nd3, contrast = c(1, 4)), "row numbers")
  expect_error(brs_marginal(m, nd3, type = "exceed", k = 11), "'k'")
})

test_that("delta-method Jacobian matches numDeriv and SEs match simulation", {
  d <- mm_test_data()
  m <- brsmm(y ~ time | time, random = ~ 1 | id, data = d, ncuts = 10)
  es <- brs_marginal(m, nd3, type = "exceed", k = 4)
  des <- .brs_marg_design(m, nd3)
  rule <- .brs_marg_rule(.brs_marg_sd_max(m, des, m$par), 8L)
  f <- function(th) rowSums(.brs_marg_pmf(th, des, m, rule)[, 5:11])
  Jn <- numDeriv::jacobian(f, m$par)
  expect_equal(unname(attr(es, "jacobian")), Jn, tolerance = 1e-7)

  set.seed(7)
  sm <- brs_marginal(m, nd3, type = "exceed", k = 4, se_method = "simulation",
                     n_sim = 4000)
  expect_equal(sm$estimate, es$estimate)
  expect_equal(sm$std.error, es$std.error, tolerance = 0.1)
  expect_true(all(sm$lower < es$estimate & sm$upper > es$estimate))
  # Logit-scale intervals stay inside (0, 1)
  expect_true(all(es$lower > 0 & es$upper < 1 & es$lower < es$estimate))
})

test_that("brs_marginal handles boundary cutpoints, sigma submodel, NA rows and levels", {
  d <- mm_test_data(G = 60, sd0 = 1.2, seed = 2)
  m <- suppressWarnings(brsmm(y ~ time | time, random = ~ 1 | id, data = d,
                              ncuts = 10, cutpoints = "boundary"))
  pm <- brs_marginal(m, nd3)
  P <- matrix(pm$estimate, nrow = 3, byrow = TRUE)
  expect_equal(rowSums(P), rep(1, 3), tolerance = 1e-12)
  # Cells of the estimated cutpoints
  cuts <- m$cut$value
  e <- .brs_marg_ends(m, m$par)
  expect_equal(e$right[1], unname(cuts[1]))
  expect_equal(e$left[11], unname(cuts[2]))
  J <- attr(pm, "jacobian")
  # Only the border scores depend on the cutpoints
  Jc <- J[, c("(cut_lo)", "(cut_hi)")]
  expect_true(all(abs(Jc[pm$score %in% c(0, 1), 1]) > 1e-6))
  expect_true(all(abs(Jc[pm$score %in% c(9, 10), 2]) > 1e-6))
  expect_true(all(abs(Jc[pm$score %in% 2:8, ]) < 1e-12))
  es <- brs_marginal(m, nd3, type = "expected_score")
  expect_true(all(is.finite(es$std.error)))

  # NA covariate row: NA output, other rows unchanged
  nd_na <- data.frame(time = c("6h", NA, "24h"))
  ena <- brs_marginal(m, nd_na, type = "expected_score")
  expect_true(is.na(ena$estimate[2]) && is.na(ena$std.error[2]))
  expect_equal(ena$estimate[c(1, 3)], es$estimate[c(1, 3)])
  # A factor with absent levels uses the levels of the fit; a new level errors
  e1 <- brs_marginal(m, data.frame(time = factor("12h")), type = "expected_score")
  expect_equal(e1$estimate, es$estimate[2])
  expect_error(brs_marginal(m, data.frame(time = "48h")), "new level")

  # brs with fixed and variable dispersion, boundary cutpoints
  f1 <- brs(y ~ time, data = d, ncuts = 10, cutpoints = "boundary")
  f2 <- brs(y ~ time | time, data = d, ncuts = 10, interval = "right")
  for (f in list(f1, f2)) {
    P <- matrix(brs_marginal(f, nd3)$estimate, nrow = 3, byrow = TRUE)
    expect_equal(rowSums(P), rep(1, 3), tolerance = 1e-12)
    expect_equal(brs_marginal(f, nd3, type = "expected_score")$estimate,
                 predict(f, newdata = nd3, type = "expected_score"), tolerance = 1e-9)
  }
  # newdata = NULL: one prediction per fitted row
  expect_equal(nrow(brs_marginal(f1, type = "expected_score")), nrow(d))
})

test_that("large random-effect SD uses an accurate rule", {
  d <- mm_test_data()
  m <- brsmm(y ~ time, random = ~ 1 | id, data = d, ncuts = 10)
  m$par[grep("re_chol", names(m$par))] <- log(4)
  P <- matrix(brs_marginal(m, nd3[1, , drop = FALSE])$estimate, nrow = 1)
  eta0 <- coef(m, "mean")[1]
  ex <- vapply(0:10, function(s) {
    stats::integrate(function(b) {
      cond_pmf(m, eta0 + b, coef(m, "precision"))[, s + 1L] * dnorm(b, 0, 4)
    }, -Inf, Inf, rel.tol = 1e-12, abs.tol = 1e-14, subdivisions = 1000L)$value
  }, numeric(1))
  expect_equal(as.vector(P), ex, tolerance = 1e-8)
})

test_that("per-group contributions add up to logLik for every integration rule", {
  d <- mm_test_data()
  for (im in c("laplace", "aghq", "qmc")) {
    m <- brsmm(y ~ time, random = ~ 1 | id, data = d, ncuts = 10, int_method = im)
    cc <- brs_loglik_contrib(m)
    expect_length(cc, 40L)
    expect_lt(abs(sum(cc) - as.numeric(logLik(m))), 1e-8)
    expect_equal(attr(cc, "n"), rep(3L, 40))
    # Held-out form on the fitting data gives the same terms
    expect_equal(brs_loglik_contrib(m, newdata = d), cc, tolerance = 1e-10)
  }
  mb <- suppressWarnings(brsmm(y ~ time | time, random = ~ 1 | id, data = d,
                               ncuts = 10, int_method = "aghq", cutpoints = "boundary"))
  expect_lt(abs(sum(brs_loglik_contrib(mb)) - as.numeric(logLik(mb))), 1e-8)
  # A different rule: AGHQ with more nodes
  c31 <- brs_loglik_contrib(mb, int_method = "aghq", n_points = 31)
  expect_lt(abs(sum(c31) - as.numeric(logLik(mb))), 1e-4)

  # brs: one term per observation
  for (f in list(brs(y ~ time, data = d, ncuts = 10),
                 brs(y ~ time | time, data = d, ncuts = 10, cutpoints = "boundary"))) {
    cf <- brs_loglik_contrib(f)
    expect_length(cf, nrow(d))
    expect_lt(abs(sum(cf) - as.numeric(logLik(f))), 1e-8)
  }
})

test_that("held-out contribution of a group matches direct integration", {
  d <- mm_test_data()
  m <- brsmm(y ~ time, random = ~ 1 | id, data = d, ncuts = 10, int_method = "aghq")
  new <- mm_test_data(G = 3, seed = 77)
  new$id <- new$id + 1000
  cc <- brs_loglik_contrib(m, newdata = new, int_method = "aghq", n_points = 31)
  expect_named(cc, as.character(1001:1003))
  sdb <- exp(coef(m, "random"))
  X <- model.matrix(~ time, new)
  eta <- drop(X %*% coef(m, "mean"))
  for (g in 1:3) {
    r <- which(new$id == 1000 + g)
    f <- function(b) {
      vapply(b, function(bb) {
        prod(cond_pmf(m, eta[r] + bb, coef(m, "precision"))[cbind(seq_along(r), new$y[r] + 1)])
      }, numeric(1)) * dnorm(b, 0, sdb)
    }
    ex <- log(stats::integrate(f, -Inf, Inf, rel.tol = 1e-12)$value)
    expect_equal(unname(cc[g]), ex, tolerance = 1e-8)
  }
})

test_that("brsmm_cv keeps groups together, matches a by-hand score and is reproducible", {
  d <- mm_test_data(G = 24)
  cv <- brsmm_cv(y ~ time, random = ~ 1 | id, data = d, k = 3, seed = 11,
                 ncuts = 10, int_method = "aghq")
  expect_s3_class(cv, "brsmm_cv")
  expect_equal(nrow(cv), 3L)
  fg <- attr(cv, "folds")[[1]]
  expect_equal(sort(names(fg)), sort(as.character(1:24)))
  expect_equal(as.vector(table(fg)), c(8L, 8L, 8L))
  expect_equal(sum(cv$n_test), nrow(d))
  expect_equal(cv$n_groups_test, c(8L, 8L, 8L))
  expect_equal(attr(cv, "log_score"), sum(cv$log_score))
  expect_equal(sum(attr(cv, "group_scores")[[1]]), sum(cv$log_score))

  # By hand: refit without each fold, integrate the held-out groups directly
  total <- 0
  for (f in 1:3) {
    test_ids <- names(fg)[fg == f]
    tr <- d[!(d$id %in% test_ids), ]
    te <- d[d$id %in% test_ids, ]
    m <- brsmm(y ~ time, random = ~ 1 | id, data = tr, ncuts = 10, int_method = "aghq")
    sdb <- exp(coef(m, "random"))
    eta <- drop(model.matrix(~ time, te) %*% coef(m, "mean"))
    for (g in test_ids) {
      r <- which(te$id == g)
      h <- function(b) {
        vapply(b, function(bb) {
          prod(cond_pmf(m, eta[r] + bb, coef(m, "precision"))[cbind(seq_along(r), te$y[r] + 1)])
        }, numeric(1)) * dnorm(b, 0, sdb)
      }
      total <- total + log(stats::integrate(h, -Inf, Inf, rel.tol = 1e-12)$value)
    }
  }
  expect_equal(attr(cv, "log_score"), total, tolerance = 1e-7)

  # Same seed, same folds and scores; the RNG state is left alone
  set.seed(5)
  before <- .Random.seed
  cv2 <- brsmm_cv(y ~ time, random = ~ 1 | id, data = d, k = 3, seed = 11,
                  ncuts = 10, int_method = "aghq")
  expect_identical(.Random.seed, before)
  expect_identical(attr(cv2, "folds"), attr(cv, "folds"))
  expect_equal(cv2$log_score, cv$log_score)
  # User folds (per row) reproduce the same assignment
  cv3 <- brsmm_cv(y ~ time, random = ~ 1 | id, data = d, folds = fg[as.character(d$id)],
                  ncuts = 10, int_method = "aghq")
  expect_equal(cv3$log_score, cv$log_score)
  expect_error(brsmm_cv(y ~ time, random = ~ 1 | id, data = d, folds = seq_len(nrow(d)),
                        ncuts = 10), "constant within each group")
})

test_that("brsmm_cv works with boundary cutpoints and records failed folds", {
  d <- mm_test_data(G = 30, sd0 = 1.2, seed = 7)
  cv <- suppressWarnings(brsmm_cv(y ~ time | time, random = ~ 1 | id, data = d, k = 3,
                                  seed = 1, ncuts = 10, cutpoints = "boundary"))
  expect_equal(nrow(cv), 3L)
  expect_true(all(is.na(cv$error)))
  expect_true(all(is.finite(cv$log_score)))
  expect_equal(attr(cv, "log_score"), sum(cv$log_score))

  # Score 10 only in one group: the fold holding it out cannot estimate c_hi
  keep <- d$id[d$y == 10][1]
  d$y[d$y == 10 & d$id != keep] <- 9
  cvf <- suppressWarnings(brsmm_cv(y ~ time, random = ~ 1 | id, data = d, k = 3,
                                   seed = 1, ncuts = 10, cutpoints = "boundary"))
  fk <- attr(cvf, "folds")[[1]][as.character(keep)]
  bad <- !is.na(cvf$error)
  expect_identical(which(bad), match(fk, cvf$fold))
  expect_match(cvf$error[bad], "not observed")
  expect_true(is.na(cvf$log_score[bad]))
  expect_true(is.na(attr(cvf, "log_score")))
  expect_true(all(is.finite(cvf$log_score[!bad])))
})

test_that("BIC(n = 'groups') counts groups; the default counts observations", {
  d <- mm_test_data()
  m <- brsmm(y ~ time, random = ~ 1 | id, data = d, ncuts = 10)
  ll <- as.numeric(logLik(m))
  expect_equal(BIC(m), -2 * ll + log(nrow(d)) * m$npar)
  expect_equal(BIC(m, n = "obs"), BIC(m))
  expect_equal(BIC(m, n = "groups"), -2 * ll + log(40) * m$npar)
})

test_that("small K and the left map", {
  set.seed(8)
  d <- expand.grid(id = 1:40, t = 1:3)
  d$y <- pmin(3, rbinom(nrow(d), 3, plogis(-0.3 + rnorm(40, sd = 0.7)[d$id])))
  m <- brsmm(y ~ 1, random = ~ 1 | id, data = d, ncuts = 3, interval = "left")
  P <- brs_marginal(m, data.frame(z = 1))$estimate
  expect_length(P, 4L)
  expect_equal(sum(P), 1, tolerance = 1e-12)
  expect_lt(abs(sum(brs_loglik_contrib(m)) - as.numeric(logLik(m))), 1e-8)
})

test_that("brsmm_cv accepts character covariates in the random-effects terms", {
  d <- mm_test_data(G = 24, sd1 = 0.4, seed = 6)
  # read.csv-style data: character columns
  d$tm <- as.character(d$time)
  d$tm2 <- ifelse(d$tm == "6h", "early", "late")
  d$time <- NULL
  cv <- suppressWarnings(brsmm_cv(y ~ tm, random = ~ 1 + tm2 | id, data = d, k = 3,
                                  seed = 2, ncuts = 10))
  expect_equal(nrow(cv), 3L)
  expect_true(all(is.na(cv$error)))
  expect_equal(attr(cv, "score_points"), 31L)
})

test_that("brsmm_cv drops incomplete rows and counts complete rows", {
  d <- mm_test_data(G = 24)
  d$y[d$id == 3] <- NA          # a group without complete rows
  d$y[d$id == 5][1] <- NA
  cv <- brsmm_cv(y ~ time, random = ~ 1 | id, data = d, k = 3, seed = 4,
                 ncuts = 10, int_method = "aghq")
  expect_equal(sum(cv$n_test), sum(!is.na(d$y)))
  expect_false("3" %in% names(attr(cv, "folds")[[1]]))
  expect_length(attr(cv, "group_scores")[[1]], 23L)
  expect_true(is.finite(attr(cv, "log_score")))
})

test_that("default held-out scoring grid stays within the node limit", {
  sp <- vapply(1:12, .brsmm_cv_score_points, integer(1))
  expect_equal(sp[1:5], c(31L, 31L, 27L, 11L, 7L))
  expect_true(all(sp >= 2L))
  expect_true(all(sp^(1:12) <= 5e5))
  d <- mm_test_data(G = 12)
  expect_error(brsmm_cv(y ~ time, random = ~ 1 + tc | id, data = d, k = 2,
                        score_points = 800, ncuts = 10), "500000-node")
})

test_that("brs_marginal warns when the cells do not partition (0, 1)", {
  d <- mm_test_data()
  m <- suppressWarnings(brsmm(y ~ time, random = ~ 1 | id, data = d, ncuts = 10,
                              lim = 0.4))
  P <- matrix(brs_marginal(m, nd3)$estimate, nrow = 3, byrow = TRUE)
  expect_true(all(rowSums(P) < 1 - 1e-3))
  expect_warning(brs_marginal(m, nd3, type = "expected_score"), "do not partition")
  expect_warning(brs_marginal(m, nd3, type = "exceed", k = 4), "do not partition")
  expect_silent(brs_marginal(brsmm(y ~ time, random = ~ 1 | id, data = d, ncuts = 10),
                             nd3, type = "expected_score"))
})

test_that("BIC with several brsmm fits returns a table with the same n rule", {
  d <- mm_test_data()
  m1 <- brsmm(y ~ time, random = ~ 1 | id, data = d, ncuts = 10)
  m2 <- brsmm(y ~ 1, random = ~ 1 | id, data = d, ncuts = 10)
  tb <- BIC(m1, m2)
  expect_s3_class(tb, "data.frame")
  expect_identical(rownames(tb), c("m1", "m2"))
  expect_equal(tb$df, c(m1$npar, m2$npar))
  expect_equal(tb$BIC, c(BIC(m1), BIC(m2)))
  tg <- BIC(m1, m2, n = "groups")
  expect_equal(tg$BIC, c(BIC(m1, n = "groups"), BIC(m2, n = "groups")))
  f <- brs(y ~ time, data = d, ncuts = 10)
  expect_equal(BIC(m1, f)$BIC, c(BIC(m1), BIC(f)))
  expect_error(BIC(m1, f, n = "groups"), "brsmm")
})

test_that("folds named by group are read by name, whatever the row order", {
  set.seed(12)
  # One row per group, rows not in level order
  d <- data.frame(id = sample(1:30), x = rnorm(30))
  d$y <- round(10 * rbeta(30, 2, 5))
  fg <- stats::setNames(rep(1:3, length.out = 30), as.character(1:30))
  ff <- .brsmm_cv_user_folds(fg, factor(d$id), by_name = TRUE)
  expect_identical(ff, fg[levels(factor(d$id))])
  # Through brsmm_cv(): the stored assignment is the named one
  d2 <- mm_test_data(G = 24)
  d2 <- d2[sample(nrow(d2)), ]
  fg2 <- stats::setNames(rep(1:3, length.out = 24), as.character(sample(1:24)))
  cv <- brsmm_cv(y ~ time, random = ~ 1 | id, data = d2, folds = fg2, ncuts = 10,
                 int_method = "aghq")
  got <- attr(cv, "folds")[[1]]
  expect_identical(got, fg2[names(got)])
  expect_equal(cv$n_test, as.vector(table(fg2[as.character(d2$id)])))
  # Neither named by the levels nor one per row
  expect_error(brsmm_cv(y ~ time, random = ~ 1 | id, data = d2, folds = 1:5,
                        ncuts = 10), "one value per row")
})

test_that("rows with a missing group are dropped as other incomplete rows", {
  d <- mm_test_data(G = 24)
  d$id[c(1, 30)] <- NA
  cv <- brsmm_cv(y ~ time, random = ~ 1 | id, data = d, k = 3, seed = 4,
                 ncuts = 10, int_method = "aghq")
  expect_equal(sum(cv$n_test), nrow(d) - 2L)
  expect_true(all(is.na(cv$error)))
})

test_that("groups without a valid random-effect mode give NA contributions", {
  d <- mm_test_data()
  m <- brsmm(y ~ time, random = ~ 1 | id, data = d, ncuts = 10)
  m$par[grep("re_chol", names(m$par))] <- -800
  expect_warning(cc <- brs_loglik_contrib(m), "no positive-definite")
  expect_true(all(is.na(cc)))
  expect_length(cc, 40L)
  expect_equal(.brs_log_penalty_cpp(), -1e6)
})

test_that("BIC tables agree whichever class comes first", {
  d <- mm_test_data()
  m <- brsmm(y ~ time, random = ~ 1 | id, data = d, ncuts = 10)
  f <- brs(y ~ time, data = d, ncuts = 10)
  a <- BIC(f, m)
  b <- BIC(m, f)
  expect_s3_class(a, "data.frame")
  expect_identical(rownames(a), c("f", "m"))
  expect_equal(a$BIC, c(BIC(f), BIC(m)))
  expect_equal(a$df, c(f$npar, m$npar))
  expect_equal(a[c("m", "f"), ], b, ignore_attr = TRUE)
  expect_equal(BIC(f), log(nrow(d)) * f$npar - 2 * f$value)
})
