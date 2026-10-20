# ============================================================================ #
# Cell-based residuals: randomized quantile residuals at given linear
# predictors (conditional, marginal or posterior random effects) and deviance
# residuals from the cell probabilities.
# ============================================================================ #

# Linear predictors of a fit: fixed part of the first parameter and the second
# parameter, one value per row.
.brs_res_eta <- function(object) {
  p <- object$p
  par <- unname(object$par)
  mm <- object$model_matrices
  eta_mu <- drop(mm$X %*% par[seq_len(p)])
  eta_phi <- if (is.null(mm$Z)) {
    rep(par[p + 1L], length(eta_mu))
  } else {
    drop(mm$Z %*% par[p + seq_len(ncol(mm$Z))])
  }
  list(mu = eta_mu, phi = eta_phi)
}

# Rows whose stored cell is the cell of an integer score of the fit's map
# (score data); NA for the other rows.
.brs_res_scores <- function(object, ends) {
  Y <- object$Y
  K <- as.integer(object$ncuts)
  s <- as.numeric(Y[, "y"])
  ok <- as.integer(Y[, "delta"]) != 0L & is.finite(s) & s == round(s) & s >= 0 & s <= K
  si <- ifelse(ok, s, 0) + 1
  ok <- ok & abs(Y[, "left"] - ends$left[si]) < 1e-8 &
    abs(Y[, "right"] - ends$right[si]) < 1e-8
  ifelse(ok, s, NA_real_)
}

# Distribution pieces of each row averaged over nodes: `eta` is n x m (one
# column per node of the random-effect integral), `w` the m weights. Returns
# the cdf F and survival S at the cell ends (exact rows: at the value), the
# log probability of the cell (log density for exact rows), E[Y] and, for
# score rows, the expected score.
.brs_res_parts <- function(object, eta, eta_phi, w = 1, need_es = TRUE) {
  Y <- object$Y
  left <- as.numeric(Y[, "left"])
  right <- as.numeric(Y[, "right"])
  yt <- as.numeric(Y[, "yt"])
  delta <- as.integer(Y[, "delta"])
  n <- length(delta)
  K <- as.integer(object$ncuts)
  eta <- matrix(eta, nrow = n)
  ends <- .brs_marg_ends(object, object$par)
  score <- .brs_res_scores(object, ends)
  # Expected scores only where needed (cell deviance signs)
  sr <- !is.na(score) & need_es
  acc <- list(Fl = 0, Fu = 0, Sl = 0, Su = 0, ey = 0, es_num = 0, es_den = 0)
  LP <- matrix(0, n, ncol(eta))
  for (j in seq_len(ncol(eta))) {
    sh <- .brs_shapes_cpp(eta[, j], eta_phi, object$link, object$link_phi,
                          object$repar)
    a <- sh$a
    b <- sh$b
    lo <- ifelse(delta == 0L, yt, left)
    hi <- ifelse(delta == 0L, yt, right)
    Fl <- stats::pbeta(lo, a, b)
    Sl <- stats::pbeta(lo, a, b, lower.tail = FALSE)
    Fu <- stats::pbeta(hi, a, b)
    Su <- stats::pbeta(hi, a, b, lower.tail = FALSE)
    # Left- and right-censored rows reach the ends of (0, 1)
    Fl[delta == 1L] <- 0
    Sl[delta == 1L] <- 1
    Fu[delta == 2L] <- 1
    Su[delta == 2L] <- 0
    acc$Fl <- acc$Fl + w[j] * Fl
    acc$Fu <- acc$Fu + w[j] * Fu
    acc$Sl <- acc$Sl + w[j] * Sl
    acc$Su <- acc$Su + w[j] * Su
    acc$ey <- acc$ey + w[j] * a / (a + b)
    LP[, j] <- .brs_obs_loglik(delta, left, right, yt, a, b)
    if (any(sr)) {
      P <- .brs_cell_probs(a[sr], b[sr], ends$left, ends$right, K)
      acc$es_num <- acc$es_num + w[j] * drop(P %*% (0:K))
      acc$es_den <- acc$es_den + w[j] * rowSums(P)
    }
  }
  # log sum_j w_j exp(lp_j), stable when the probabilities are tiny
  mx <- apply(LP, 1L, max)
  lp <- if (ncol(LP) == 1L) LP[, 1L] else mx + log(drop(exp(LP - mx) %*% w))
  es <- rep(NA_real_, n)
  if (any(sr)) es[sr] <- acc$es_num / acc$es_den
  list(Fl = acc$Fl, Fu = acc$Fu, Sl = acc$Sl, Su = acc$Su, lp = lp,
       ey = acc$ey, score = score, es = es, delta = delta, yt = yt)
}

# Randomized quantile residuals from the cell cdf/survival pieces: uniform on
# (F(l), F(u)), drawn on the survival scale in the upper half so that upper
# tail cells keep their precision; one runif() per row.
.brs_rqr_parts <- function(parts) {
  n <- length(parts$Fl)
  u <- stats::runif(n)
  upper <- parts$Fl > 0.5
  v <- parts$Fl + u * pmax(parts$Fu - parts$Fl, 0)
  sv <- parts$Sl - u * pmax(parts$Sl - parts$Su, 0)
  v <- pmin(pmax(v, 1e-10), 1 - 1e-10)
  sv <- pmin(pmax(sv, 1e-10), 1 - 1e-10)
  ifelse(upper, stats::qnorm(sv, lower.tail = FALSE), stats::qnorm(v))
}

# Deviance residuals of censored rows from the cell probabilities:
# sign * sqrt(-2 log P), with sign +1 when the observed score is at least the
# expected score (cell center vs E[Y] for non-score rows); NA for exact rows.
.brs_cell_dev_parts <- function(parts) {
  diff <- ifelse(is.na(parts$score), parts$yt - parts$ey, parts$score - parts$es)
  sgn <- ifelse(diff >= 0, 1, -1)
  d <- sgn * sqrt(pmax(-2 * parts$lp, 0))
  d[parts$delta == 0L] <- NA_real_
  d
}

# Nodes of the marginal (population-averaged) distribution of each row of a
# brsmm fit: eta + sd_u z on the rule of brs_marginal(), u = x_r' b.
.brsmm_res_marginal_nodes <- function(object, eta_mu, n_points = 8L) {
  theta <- object$par
  Xr <- object$model_matrices$Xr
  sd_u <- sqrt(rowSums((Xr %*% object$random$L)^2))
  gh <- .brs_marg_rule(max(sd_u), n_points)
  list(eta = eta_mu + outer(sd_u, gh$x), w = gh$w)
}

# Log posterior (up to a constant) of the random effects of every group, as a
# function of a G x q matrix B (one row per group).
.brsmm_logpost_fun <- function(object, eta_mu, eta_phi) {
  Y <- object$Y
  left <- as.numeric(Y[, "left"])
  right <- as.numeric(Y[, "right"])
  yt <- as.numeric(Y[, "yt"])
  delta <- as.integer(Y[, "delta"])
  Xr <- object$model_matrices$Xr
  gi <- object$group_index
  Linv <- solve(object$random$L)
  function(B) {
    eta <- eta_mu + rowSums(Xr * B[gi, , drop = FALSE])
    sh <- .brs_shapes_cpp(eta, eta_phi, object$link, object$link_phi, object$repar)
    lp <- .brs_obs_loglik(delta, left, right, yt, sh$a, sh$b)
    as.numeric(rowsum(lp, gi, reorder = TRUE)) - 0.5 * rowSums((B %*% t(Linv))^2)
  }
}

# Negative Hessian of `logpost` at the rows of B0 by central differences
# (step h), all groups at once: a G x q x q array.
.brsmm_curvature <- function(logpost, B0, h = 1e-3) {
  G <- nrow(B0)
  q <- ncol(B0)
  H <- array(0, c(G, q, q))
  f0 <- logpost(B0)
  at <- function(k, sk, l = NULL, sl = 0) {
    B <- B0
    B[, k] <- B[, k] + sk
    if (!is.null(l)) B[, l] <- B[, l] + sl
    logpost(B)
  }
  for (k in seq_len(q)) {
    H[, k, k] <- -(at(k, h) - 2 * f0 + at(k, -h)) / h^2
    if (k < q) {
      for (l in (k + 1L):q) {
        H[, k, l] <- H[, l, k] <- -(at(k, h, l, h) - at(k, h, l, -h) -
          at(k, -h, l, h) + at(k, -h, l, -h)) / (4 * h^2)
      }
    }
  }
  H
}

# One draw of the random effects of each group from its approximate
# posterior: sampling-importance-resampling with `n_draws` proposals from
# N(mode, H^-1) (H: curvature of the log posterior at the mode), weights
# posterior / proposal. Returns a G x q matrix.
.brsmm_posterior_b <- function(object, eta_mu, eta_phi, n_draws) {
  G <- length(object$random$levels)
  L <- object$random$L
  q <- ncol(L)
  B0 <- matrix(object$random$mode_b, nrow = G)
  logpost <- .brsmm_logpost_fun(object, eta_mu, eta_phi)
  H <- .brsmm_curvature(logpost, B0)
  # Proposal factor C_g (C C' = H_g^-1); the prior factor L when H_g is not
  # positive definite
  C <- array(0, c(G, q, q))
  for (g in seq_len(G)) {
    R <- tryCatch(chol(matrix(H[g, , ], q, q)), error = function(e) NULL)
    C[g, , ] <- if (is.null(R) || any(!is.finite(R))) L else backsolve(R, diag(q))
  }
  Zs <- array(stats::rnorm(n_draws * G * q), c(G, q, n_draws))
  logw <- matrix(0, G, n_draws)
  props <- vector("list", n_draws)
  for (r in seq_len(n_draws)) {
    z <- matrix(Zs[, , r], G, q)
    B <- B0
    for (k in seq_len(q)) {
      B[, k] <- B0[, k] + rowSums(matrix(C[, k, ], G, q) * z)
    }
    props[[r]] <- B
    # Proposal log density up to a per-group constant: -|z|^2 / 2
    logw[, r] <- logpost(B) + 0.5 * rowSums(z^2)
  }
  wt <- exp(logw - apply(logw, 1L, max))
  wt <- wt / rowSums(wt)
  u <- stats::runif(G)
  pick <- pmin(rowSums(u > t(apply(wt, 1L, cumsum))) + 1L, n_draws)
  matrix(vapply(seq_len(G), function(g) props[[pick[g]]][g, ], numeric(q)),
         nrow = G, ncol = q, byrow = TRUE)
}

# Residuals of type "rqr" or "cell_deviance" of a brsmm fit for the chosen
# treatment of the random effects.
.brsmm_res_cell <- function(object, type, re, n_draws) {
  need_es <- type == "cell_deviance"
  et <- .brs_res_eta(object)
  Xr <- object$model_matrices$Xr
  G <- length(object$random$levels)
  if (re == "marginal") {
    nd <- .brsmm_res_marginal_nodes(object, et$mu)
    parts <- .brs_res_parts(object, nd$eta, et$phi, nd$w, need_es = need_es)
  } else {
    B <- if (re == "posterior") {
      .brsmm_posterior_b(object, et$mu, et$phi, n_draws)
    } else {
      matrix(object$random$mode_b, nrow = G)
    }
    eta <- et$mu + rowSums(Xr * B[object$group_index, , drop = FALSE])
    parts <- .brs_res_parts(object, eta, et$phi, need_es = need_es)
  }
  if (type == "rqr") .brs_rqr_parts(parts) else .brs_cell_dev_parts(parts)
}
