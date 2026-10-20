# Estimated boundary cutpoints (cutpoints = "boundary")

# NRS-11-like scores: 80 patients x 3 times, random intercept, extra covariate z
cut_test_data <- function() {
  set.seed(2023)
  nrs <- expand.grid(id = 1:80, time = c("6h", "12h", "24h"))
  nrs$group <- factor(paste0("g", (nrs$id - 1) %% 4 + 1))
  nrs$z <- rnorm(nrow(nrs))
  b0 <- rnorm(80, sd = 0.6)
  eta <- -1.3 + b0[nrs$id] + c(0, 0.75, 0.3)[nrs$time] +
    c(0, -0.1, 0.05, 0.1)[nrs$group]
  shp <- brs_repar(mu = plogis(eta), phi = 0.3, repar = 2)
  nrs$y <- round(10 * rbeta(nrow(nrs), shp$shape1, shp$shape2))
  nrs
}

# Scores coarsened with known boundary cutpoints (K = 10 unless given)
cut_sim <- function(n, cuts = c(0.11, 0.9), K = 10L, interval = "mid", seed = 1) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n))
  brs_sim(~ x, data = d, beta = c(-0.4, 0.6), phi = qlogis(0.25), ncuts = K,
          interval = interval, cutpoints = cuts)
}

test_that("cutpoints = 'fixed' reproduces the fixed-map fits", {
  nrs <- cut_test_data()
  f1 <- brs(y ~ time + group, data = nrs, ncuts = 10)
  f2 <- brs(y ~ time | z, data = nrs, ncuts = 10, interval = "right")
  f3 <- brs(y ~ time, data = nrs, ncuts = 10, repar = 1)
  m1 <- suppressWarnings(brsmm(y ~ time, random = ~ 1 | id, data = nrs, ncuts = 10))
  ref <- list(
    f1 = list(value = -507.742212721466, par = c(
      -1.21754242846245, 0.612065481342978, 0.319160325962717,
      0.067372255573696, -0.23408541944632, 0.0740216841704957, -0.576939455164775
    ), se = c(0.208696820354, 0.199187637744, 0.202997304657, 0.229467222065,
              0.231091778996, 0.227617657001, 0.0961381172548)),
    f2 = list(value = -510.824850824695, par = c(
      -1.13283014913746, 0.556975666680916, 0.305318032613135,
      -0.839824179884508, -0.000472843215514938
    ), se = c(0.139623929615, 0.186538790441, 0.189731147714, 0.0980977198767,
              0.100186913256)),
    f3 = list(value = -508.949681912759, par = c(
      -1.23804550855672, 0.604606994936061, 0.323911202815055, 0.562534346392089
    ), se = c(0.150493264751, 0.199937424866, 0.203133066946, 0.0958247023985)),
    m1 = list(value = -508.245050083974, par = c(
      -1.27638467887282, 0.625904295906637, 0.340596003858071,
      -0.634022533308435, -1.06788317196871
    ), se = c(0.158371412496, 0.199584757948, 0.202445329016, 0.113041475903,
              0.470568870445))
  )
  fits <- list(f1 = f1, f2 = f2, f3 = f3, m1 = m1)
  for (nm in names(fits)) {
    f <- fits[[nm]]
    expect_equal(f$value, ref[[nm]]$value, tolerance = 1e-8)
    expect_equal(unname(f$par), ref[[nm]]$par, tolerance = 1e-5)
    expect_equal(unname(sqrt(diag(vcov(f)))), ref[[nm]]$se, tolerance = 1e-5)
    expect_identical(f$cutpoints, "fixed")
    expect_null(f$cut)
  }
  # The explicit default is the same computation
  f1b <- brs(y ~ time + group, data = nrs, ncuts = 10, cutpoints = "fixed")
  expect_identical(f1b$par, f1$par)
  expect_identical(f1b$hessian, f1$hessian)
  expect_identical(f1b$value, f1$value)
})

test_that("analytic cutpoint gradient matches numDeriv for all maps", {
  nrs <- cut_test_data()
  for (iv in c("mid", "right", "left")) {
    for (fm in list(y ~ time + group, y ~ time | z)) {
      fb <- brs(fm, data = nrs, ncuts = 10, interval = iv, cutpoints = "boundary")
      ob <- .brs_cut_objective(fb)
      set.seed(7)
      par <- fb$par + rnorm(length(fb$par), sd = 0.15)
      ga <- ob$gr(par)
      gn <- numDeriv::grad(ob$fn, par)
      k <- length(par) - 1:0
      expect_lt(max(abs(ga[k] - gn[k]) / pmax(abs(gn[k]), 1)), 1e-6)
      expect_lt(max(abs(ga - gn) / pmax(abs(gn), 1)), 1e-5)
    }
  }
})

test_that("boundary Hessian and standard errors match numDeriv", {
  nrs <- cut_test_data()
  for (fm in list(y ~ time + group, y ~ time | z)) {
    fb <- brs(fm, data = nrs, ncuts = 10, cutpoints = "boundary")
    ob <- .brs_cut_objective(fb)
    Hn <- numDeriv::hessian(ob$fn, fb$par)
    expect_lt(max(abs(fb$hessian - Hn)) / max(abs(Hn)), 1e-5)
    se_n <- sqrt(diag(solve(-Hn)))
    expect_equal(unname(sqrt(diag(vcov(fb)))), se_n, tolerance = 1e-5)
    # hessian_method = "numDeriv" uses the same objective
    fn <- brs(fm, data = nrs, ncuts = 10, cutpoints = "boundary",
              hessian_method = "numDeriv")
    expect_equal(unname(sqrt(diag(vcov(fn)))), se_n, tolerance = 1e-4)
  }
})

test_that("score probabilities sum to one under boundary cutpoints", {
  for (iv in c("mid", "right", "left")) {
    for (K in c(2L, 3L, 4L, 10L)) {
      s <- cut_sim(600, cuts = NULL, K = K, interval = iv, seed = K)
      fb <- suppressWarnings(brs(y ~ x, data = s, cutpoints = "boundary"))
      P <- brs_predict_scoreprob(fb)
      expect_equal(unname(rowSums(P)), rep(1, nrow(P)), tolerance = 1e-8)
      expect_true(all(P >= 0))
      # Expected score uses the same probabilities
      expect_equal(predict(fb, type = "expected_score"),
                   as.numeric(P %*% (0:K)), tolerance = 1e-12)
      cv <- fb$cut$value
      expect_true(cv[1] > 0 && cv[1] < cv[2] && cv[2] < 1)
    }
  }
})

test_that("the fixed map is nested in the boundary model", {
  nrs <- cut_test_data()
  for (iv in c("mid", "right")) {
    for (fm in list(y ~ time + group, y ~ time | z)) {
      f0 <- brs(fm, data = nrs, ncuts = 10, interval = iv)
      fb <- brs(fm, data = nrs, ncuts = 10, interval = iv, cutpoints = "boundary")
      expect_gte(fb$value, f0$value - 1e-8)
      # Boundary likelihood at the map's g reproduces the fixed likelihood
      ob <- .brs_cut_objective(fb)
      expect_equal(ob$fn(c(f0$par, fb$cut$start)), f0$value, tolerance = 1e-10)
      expect_equal(.brs_cut_values(fb$cut$start, 10L, .brs_cut_anchors(10L, iv))$value,
                   fb$cut$fixed, tolerance = 1e-14)
      a <- anova(f0, fb)
      expect_equal(a$`Chi Df`[2], 2)
      expect_equal(a$Chisq[2], 2 * (fb$value - f0$value), tolerance = 1e-10)
      expect_true(any(grepl("boundary", attr(a, "heading"))))
    }
  }
  # Map values of g: logit(1/3), logit(2/3) under mid; 0, 0 under right/left
  expect_equal(.brs_cut_par(c(0.05, 0.95), 10L, .brs_cut_anchors(10L, "mid")),
               qlogis(c(1 / 3, 2 / 3)))
  expect_equal(.brs_cut_par(c(1, 10) / 11, 10L, .brs_cut_anchors(10L, "right")),
               c(0, 0))
})

test_that("boundary cutpoints are recovered from simulated data", {
  s <- cut_sim(2000, cuts = c(0.11, 0.9))
  fb <- brs(y ~ x, data = s, cutpoints = "boundary")
  tab <- summary(fb)$cutpoints
  z <- (tab[, "cutpoint"] - c(0.11, 0.9)) / tab[, "SE(cut)"]
  expect_true(all(abs(z) < 4))
  z_beta <- (coef(fb, model = "mean") - c(-0.4, 0.6)) / sqrt(diag(vcov(fb, model = "mean")))
  expect_true(all(abs(z_beta) < 4))
  # Score scale is K times the unit scale under mid
  expect_equal(tab[, "score"], 10 * tab[, "cutpoint"])
  expect_equal(unname(coef(fb, model = "cutpoints")), unname(fb$cut$par))
  expect_equal(dim(vcov(fb, model = "cutpoints")), c(2L, 2L))
  expect_equal(nrow(confint(fb, model = "cutpoints")), 2L)
  expect_equal(attr(logLik(fb), "df"), 5)
  expect_equal(AIC(fb), -2 * fb$value + 10)
  expect_equal(nobs(fb), 2000)
  out <- capture.output(print(fb), summary(fb))
  expect_true(any(grepl("Boundary cutpoints", out)))
})

test_that("small scales K = 2, 3, 4 fit with ordered cutpoints", {
  for (K in 2:4) {
    s <- cut_sim(1500, cuts = NULL, K = K, seed = 10 + K)
    fb <- suppressWarnings(brs(y ~ x, data = s, cutpoints = "boundary"))
    cv <- fb$cut$value
    expect_true(cv[1] > 0 && cv[1] < cv[2] && cv[2] < 1)
    if (K >= 3L) {
      a <- .brs_cut_anchors(K, "mid")
      expect_lt(cv[1], a$u1)
      expect_gt(cv[2], a$l_km1)
    }
    ob <- .brs_cut_objective(fb)
    gn <- numDeriv::grad(ob$fn, fb$par + 0.05)
    ga <- ob$gr(fb$par + 0.05)
    k <- length(gn) - 1:0
    expect_lt(max(abs(ga[k] - gn[k]) / pmax(abs(gn[k]), 1)), 1e-6)
    # Stored cells: score 1 starts at c_lo, score K - 1 ends at c_hi
    Y <- fb$Y
    expect_equal(unique(Y[Y[, "y"] == 1, "left"]), cv[[1]])
    expect_equal(unique(Y[Y[, "y"] == K - 1, "right"]), cv[[2]])
    expect_equal(unique(Y[Y[, "y"] == 0, "right"]), cv[[1]])
    expect_equal(unique(Y[Y[, "y"] == K, "left"]), cv[[2]])
  }
  # K = 2: the stick-breaking order and its Jacobian
  g <- c(0.3, -0.4)
  v <- .brs_cut_values(g, 2L, .brs_cut_anchors(2L, "mid"))
  num <- numDeriv::jacobian(function(gg) .brs_cut_values(gg, 2L, NULL)$value, g)
  expect_equal(v$jacobian, num, tolerance = 1e-8)
  expect_equal(unname(.brs_cut_par(v$value, 2L, NULL)), g)
})

test_that("boundary cutpoints stop on data they cannot use", {
  base <- data.frame(y = c(0, 1, 9, 10, 4, 5, 0, 10), x = 1:8)
  # No score K: c_hi not identified
  expect_error(brs(y ~ x, data = transform(base, y = pmin(y, 9)), ncuts = 10,
                   cutpoints = "boundary"), "not observed")
  # Exact observations
  expect_error(suppressWarnings(brs(y ~ x, data = rbind(base, data.frame(y = 0.5, x = 9)),
                                    ncuts = 10, cutpoints = "boundary")), "exact")
  # lim != 0.5
  expect_error(suppressWarnings(brs(y ~ x, data = base, ncuts = 10, lim = 0.4,
                                    cutpoints = "boundary")), "lim = 0.5")
  # Analyst-supplied censoring (brs_prep Mode 2)
  pp <- suppressWarnings(suppressMessages(brs_prep(
    transform(base, delta = c(NA, NA, NA, NA, 1, NA, NA, NA)), ncuts = 10)))
  expect_error(brs(y ~ x, data = pp, cutpoints = "boundary"), "analyst")
  # Analyst bounds (Mode 3)
  pp3 <- suppressWarnings(suppressMessages(brs_prep(data.frame(
    y = c(0, 1, 9, 10, NA), left = c(NA, NA, NA, NA, 2), right = c(NA, NA, NA, NA, 5),
    x = 1:5), ncuts = 10)))
  expect_error(brs(y ~ x, data = pp3, cutpoints = "boundary"))
  expect_error(brs(y ~ x, data = base, ncuts = 10, cutpoints = "bnd"), "fixed")
  # coef(model = "cutpoints") of a fixed-map fit
  f0 <- suppressWarnings(brs(y ~ 1, data = base, ncuts = 10))
  expect_error(coef(f0, model = "cutpoints"), "no estimated cutpoints")
  # A collapsing border cell is reported
  cs <- .brs_cut_setup(brs_check(base$y, ncuts = 10), 10L, 0.5, "mid")
  expect_warning(.brs_cut_component(cs, c(12, 0)), "edge of their range")
  expect_silent(.brs_cut_component(cs, cs$start))
  # brs_sim() checks the cutpoints against the map
  d <- data.frame(x = 1:10)
  expect_error(brs_sim(~ x, data = d, beta = c(0, 0.1), ncuts = 10,
                       cutpoints = c(0.2, 0.9)), "c_lo")
  expect_error(brs_sim(~ x, data = d, beta = c(0, 0.1), ncuts = 10,
                       cutpoints = c(0.1, 0.9), delta = 3L), "delta = NULL")
})

test_that("brs_sim() coarsens with the given boundary cutpoints", {
  set.seed(5)
  d <- data.frame(x = rnorm(3000))
  set.seed(6)
  s <- brs_sim(~ x, data = d, beta = c(0, 0.5), phi = qlogis(0.3), ncuts = 10,
               cutpoints = c(0.11, 0.9))
  set.seed(6)
  s0 <- brs_sim(~ x, data = d, beta = c(0, 0.5), phi = qlogis(0.3), ncuts = 10)
  # Same draws: only rows in the moved border cells change score
  expect_true(all(abs(s$y - s0$y) <= 1))
  expect_true(all(s$y[s0$y %in% 2:8] == s0$y[s0$y %in% 2:8]))
  expect_gt(sum(s$y == 0), sum(s0$y == 0))
  expect_equal(.brs_score_from_unit_cut(c(0.1, 0.11, 0.149, 0.899, 0.9), 10L, "mid",
                                        c(0.11, 0.9)), c(0, 1, 1, 9, 10))
})

test_that("bootstrap and cross-validation use the estimated cutpoints", {
  s <- cut_sim(600, cuts = c(0.12, 0.88), seed = 3)
  fb <- brs(y ~ x, data = s, cutpoints = "boundary")
  set.seed(1)
  bt <- suppressWarnings(brs_bootstrap(fb, R = 10L))
  expect_equal(bt$parameter, names(fb$par))
  expect_true(all(is.finite(bt$se_boot)))
  setup <- .brs_boot_setup(fb)
  expect_equal(setup$cuts, fb$cut$value)
  expect_true(all(setup$mech == "score"))
  set.seed(2)
  cvr <- brs_cv(y ~ x, data = s, k = 3, cutpoints = "boundary")
  expect_true(all(cvr$converged))
  # Held-out log score equals the fit's own contributions on the training data
  m <- .brs_cv_metrics(fb, s)
  expect_equal(m$log_score * nrow(s), fb$value, tolerance = 1e-8)
  # Randomized quantile residuals use the stored (estimated) cells
  set.seed(3)
  r <- residuals(fb, type = "rqr")
  expect_true(all(is.finite(r)))
})

test_that("brsmm estimates boundary cutpoints with a random intercept", {
  set.seed(42)
  G <- 60L
  d <- expand.grid(id = seq_len(G), t = 1:3)
  b <- rnorm(G, sd = 0.7)
  eta <- -0.6 + 0.3 * (d$t - 2) + b[d$id]
  shp <- brs_repar(mu = plogis(eta), phi = 0.25, repar = 2)
  ystar <- rbeta(nrow(d), shp$shape1, shp$shape2)
  d$y <- .brs_score_from_unit_cut(ystar, 10L, "mid", c(0.11, 0.9))
  m0 <- brsmm(y ~ t, random = ~ 1 | id, data = d, ncuts = 10, int_method = "aghq",
              n_points = 9)
  for (im in c("laplace", "aghq", "qmc")) {
    mb <- brsmm(y ~ t, random = ~ 1 | id, data = d, ncuts = 10, int_method = im,
                n_points = 9, qmc_points = 256, cutpoints = "boundary")
    expect_equal(names(mb$par)[length(mb$par) - 1:0], c("(cut_lo)", "(cut_hi)"))
    expect_true(mb$diagnostics$hessian_nd)
    expect_true(all(is.finite(sqrt(diag(vcov(mb))))))
    expect_equal(attr(logLik(mb), "df"), 6)
    tab <- summary(mb)$cutpoints
    expect_equal(dim(tab), c(2L, 7L))
    expect_true(all(abs(tab[, "cutpoint"] - c(0.11, 0.9)) < 4 * tab[, "SE(cut)"] + 0.02))
    P <- .brs_score_prob_matrix(mb$fitted_mu, mb$fitted_phi, 2L, 10L, 0.5, 0:10,
                                "mid", cuts = .brs_cuts_of(mb))
    expect_equal(unname(rowSums(P)), rep(1, nrow(P)), tolerance = 1e-8)
    expect_equal(predict(mb, type = "expected_score"), as.numeric(P %*% (0:10)),
                 tolerance = 1e-10)
    if (im == "aghq") {
      expect_gte(mb$value, m0$value - 1e-6)
      # The Hessian matches numDeriv on the marginal likelihood
      fn <- function(par) {
        cs <- .brs_cut_setup(.brs_cut_fixed_Y(mb), 10L, 0.5, "mid")
        e <- .brs_cut_endpoints(.brs_cut_values(par[5:6], 10L, cs$anchors), cs)
        .brsmm_loglik_eigen(par[1:4], mb$model_matrices$X, mb$model_matrices$Z,
                            mb$model_matrices$Xr, e$left, e$right,
                            as.numeric(mb$Y[, "yt"]), mb$delta, mb$group_index,
                            0L, 0L, 2L, 1L, 9L)
      }
      expect_equal(fn(mb$par), mb$value, tolerance = 1e-8)
      Hn <- numDeriv::hessian(fn, unname(mb$par))
      expect_equal(unname(sqrt(diag(solve(-mb$hessian)))), sqrt(diag(solve(-Hn))),
                   tolerance = 1e-3)
      a <- anova(m0, mb)
      expect_equal(a$`Chi Df`[2], 2)
    }
  }
})
