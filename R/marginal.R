# ============================================================================ #
# Marginal (population-averaged) score predictions
# ============================================================================ #

#' Population-averaged score probabilities, expected scores and exceedance
#'
#' @description
#' Score distribution of a \code{\link{brs}} or \code{\link{brsmm}} fit at
#' given covariate values, with standard errors and intervals: the
#' probabilities \eqn{P(S = s)} of the scores \eqn{0, \ldots, K}, the
#' expected score or the exceedance probabilities \eqn{P(S \ge k)}, and
#' differences of these between pairs of rows. For \code{brsmm} fits the
#' random effects are integrated out (marginal, population-averaged
#' quantities) or set to zero (\code{conditional = TRUE}).
#'
#' @details
#' For a row with fixed-effect linear predictor \eqn{\eta} and random-effect
#' design \eqn{x_r}, the marginal score probabilities are
#' \deqn{P(S = s) = \int P(S = s \mid \eta + x_r^\top b)\,
#'   \varphi(b; 0, D)\, db = E\{P(S = s \mid \eta + u)\}, \qquad
#'   u \sim N(0, x_r^\top D x_r),}
#' because \eqn{b} enters the model only through \eqn{u = x_r^\top b} (the
#' dispersion submodel has no random effects), so the integral is
#' one-dimensional whatever the number of random effects. It is computed on
#' \eqn{z = u / \mathrm{SD}(u) \in [-8.5, 8.5]} by a composite
#' Gauss-Legendre rule with \code{n_points} nodes per panel and panels of
#' width at most 1 on the scale of \eqn{z} and at most 0.5 on the scale of
#' \eqn{u} (the conditional probabilities can change quickly with \eqn{u}
#' when the precision is high or the variance large, where Gauss-Hermite
#' rules converge slowly). A warning is given when halving the panel width
#' changes a value by more than \eqn{10^{-6}}.
#'
#' \eqn{P(S = s \mid \cdot)} is the probability of the cell of \eqn{s} under
#' the latent beta law, with the cells of the fit's \code{interval} map or its
#' estimated boundary cutpoints (as in \code{\link{brs_predict_scoreprob}},
#' but without a probability floor), so each row sums to 1 whenever the cells
#' partition (0, 1) (\code{lim = 0.5}, or \code{interval = "right"} or
#' \code{"left"}). \code{"expected_score"} is \eqn{\sum_s s\, P(S = s)} on the
#' \eqn{0, \ldots, K} scale and \code{"exceed"} is
#' \eqn{P(S \ge k) = \sum_{s \ge k} P(S = s)}.
#'
#' Under a nonlinear link the marginal quantities differ from the
#' conditional ones at \eqn{b = 0}: the latter describe a subject whose
#' random effects are zero, not the population average, and the gap grows
#' with the random-effect variance. For a \code{brs} fit both coincide.
#'
#' Uncertainty: with \code{se_method = "delta"} (default), delta-method
#' standard errors from a Richardson central-difference Jacobian of the
#' quantities in all parameters of \code{coef(object)} (regression,
#' precision, random-effect Cholesky and boundary-cutpoint parameters) and
#' \code{vcov(object)}. The Wald intervals are computed on the logit scale for
#' probabilities and for \eqn{E[S]/K} and transformed back (so they stay in
#' range), and on the natural scale for contrasts. With
#' \code{se_method = "simulation"}, \code{n_sim} parameter vectors are drawn
#' from \eqn{N(\hat\theta, V)}; the standard error is the standard deviation
#' of the quantity over the draws and the interval its percentile interval
#' (the draws use the RNG: call \code{set.seed()} first for reproducible
#' values). The estimate is always the plug-in value at \eqn{\hat\theta}.
#'
#' Rows of \code{newdata} with missing covariates give \code{NA}; the response
#' and the grouping variable are not needed.
#'
#' @param object A fitted \code{"brs"} or \code{"brsmm"} object.
#' @param newdata Data frame of covariate values, one prediction per row;
#'   \code{NULL} (default) uses the rows of the fit.
#' @param type \code{"prob"} (default; probabilities of the scores
#'   \eqn{0, \ldots, K}), \code{"expected_score"} or \code{"exceed"}.
#' @param k Integer thresholds in \eqn{0, \ldots, K} for
#'   \code{type = "exceed"}; \code{NULL} (default) uses \eqn{1, \ldots, K}.
#' @param contrast Optional pairs of rows of \code{newdata} to compare: a
#'   two-column matrix with one pair per row, a list of length-2 vectors, or
#'   a single length-2 vector. The pair \code{c(i, j)} gives the quantity of
#'   row \code{i} minus that of row \code{j}.
#' @param conditional Logical; \code{TRUE} sets the random effects to zero
#'   instead of integrating them out (no effect for \code{brs} fits).
#' @param level Confidence level of the intervals (default 0.95).
#' @param se_method \code{"delta"} (default) or \code{"simulation"}.
#' @param n_sim Number of parameter draws for \code{se_method = "simulation"}
#'   (default 1000).
#' @param n_points Gauss-Legendre nodes per panel of the random-effect
#'   integral (default 8; see Details).
#'
#' @return A data frame of class \code{"brs_marginal"} with the columns of
#'   \code{newdata} (or, for contrasts, \code{contrast}, \code{row1} and
#'   \code{row2}), \code{row}, \code{score} (\code{type = "prob"}) or
#'   \code{k} (\code{type = "exceed"}), \code{estimate}, \code{std.error},
#'   \code{lower} and \code{upper}. Attributes: \code{type},
#'   \code{conditional}, \code{level}, \code{se_method}, \code{n_nodes}
#'   (total quadrature nodes; \code{NA} without random effects)
#'   and, for the delta method, \code{jacobian} (quantities by parameters).
#'
#' @seealso \code{\link{brsmm}}, \code{\link{brs_predict_scoreprob}},
#'   \code{\link{predict.brsmm}}, \code{\link{brs_loglik_contrib}}
#'
#' @examples
#' # Synthetic NRS-11 scores of 40 patients at three times. Simulated data.
#' set.seed(21)
#' nrs <- expand.grid(id = 1:40, time = c("6h", "12h", "24h"))
#' b0 <- rnorm(40, sd = 0.8)
#' eta <- -1.3 + c(0, 0.75, 0.3)[nrs$time] + b0[nrs$id]
#' shp <- brs_repar(mu = plogis(eta), phi = 0.2, repar = 2)
#' nrs$y <- round(10 * rbeta(nrow(nrs), shp$shape1, shp$shape2))
#' m <- brsmm(y ~ time, random = ~ 1 | id, data = nrs, ncuts = 10)
#'
#' nd <- data.frame(time = c("6h", "12h", "24h"))
#' # Population-averaged mean score and P(S >= 4) at each time
#' brs_marginal(m, nd, type = "expected_score")
#' brs_marginal(m, nd, type = "exceed", k = 4)
#' # Same at b = 0 (a typical patient): not the population average
#' brs_marginal(m, nd, type = "expected_score", conditional = TRUE)
#' # Change from 6h: 12h - 6h and 24h - 6h
#' brs_marginal(m, nd, type = "expected_score",
#'              contrast = rbind(c(2, 1), c(3, 1)))
#'
#' @export
brs_marginal <- function(object,
                         newdata = NULL,
                         type = c("prob", "expected_score", "exceed"),
                         k = NULL,
                         contrast = NULL,
                         conditional = FALSE,
                         level = 0.95,
                         se_method = c("delta", "simulation"),
                         n_sim = 1000L,
                         n_points = 8L) {
  if (!inherits(object, c("brs", "brsmm"))) {
    stop("'object' must be a \"brs\" or \"brsmm\" fit.", call. = FALSE)
  }
  type <- match.arg(type)
  se_method <- match.arg(se_method)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("'level' must be a number in (0, 1).", call. = FALSE)
  }
  if (!is.logical(conditional) || length(conditional) != 1L || is.na(conditional)) {
    stop("'conditional' must be TRUE or FALSE.", call. = FALSE)
  }
  n_points <- as.integer(n_points)
  if (length(n_points) != 1L || is.na(n_points) || n_points < 2L || n_points > 50L) {
    stop("'n_points' must be an integer in [2, 50].", call. = FALSE)
  }
  n_sim <- as.integer(n_sim)
  if (se_method == "simulation" && (length(n_sim) != 1L || is.na(n_sim) || n_sim < 10L)) {
    stop("'n_sim' must be an integer >= 10.", call. = FALSE)
  }
  K <- as.integer(object$ncuts)
  if (type == "exceed") {
    if (is.null(k)) k <- seq_len(K)
    if (!is.numeric(k) || !length(k) || anyNA(k) || any(k != round(k)) ||
        any(k < 0) || any(k > K)) {
      stop("'k' must be integers in [0, ", K, "].", call. = FALSE)
    }
    k <- as.integer(k)
  }

  des <- .brs_marg_design(object, newdata)
  n <- des$n
  # Distinct covariate rows only: the pmf is computed once per pattern
  key <- do.call(paste, c(as.data.frame(cbind(des$X, des$Z, des$Xr)), sep = "\r"))
  first <- !duplicated(key)
  map <- match(key, key[first])
  des_u <- des
  des_u$X <- des$X[first, , drop = FALSE]
  if (!is.null(des$Z)) des_u$Z <- des$Z[first, , drop = FALSE]
  if (!is.null(des$Xr)) des_u$Xr <- des$Xr[first, , drop = FALSE]
  pairs <- .brs_marg_pairs(contrast, n)
  mixed <- inherits(object, "brsmm") && !conditional
  theta <- object$par
  # Quadrature rule fixed at the estimate (kept while differentiating)
  sd_max <- if (mixed) .brs_marg_sd_max(object, des_u, theta) else 0
  gh <- .brs_marg_rule(sd_max, n_points)

  n_col <- switch(type, prob = K + 1L, expected_score = 1L, exceed = length(k))
  # Quantities of the complete rows (one row per newdata row, one column per
  # score or threshold), then the requested contrasts
  quant <- function(theta, gh_rule = gh) {
    P <- .brs_marg_pmf(theta, des_u, object, gh_rule)[map, , drop = FALSE]
    Q <- switch(type,
      prob = P,
      expected_score = P %*% (0:K),
      exceed = vapply(k, function(kk) {
        rowSums(P[, (kk:K) + 1L, drop = FALSE])
      }, numeric(nrow(P)))
    )
    Q <- matrix(Q, nrow = nrow(P), ncol = n_col)
    full <- matrix(NA_real_, n, ncol(Q))
    full[des$ok, ] <- Q
    if (!is.null(pairs)) {
      full <- full[pairs[, 1L], , drop = FALSE] - full[pairs[, 2L], , drop = FALSE]
    }
    full
  }

  est <- quant(theta)
  if (sd_max > 0 && any(des$ok)) {
    # Quadrature check against a finer rule
    alt <- quant(theta, .brs_marg_rule(sd_max, n_points, refine = TRUE))
    dev <- max(abs(alt - est), na.rm = TRUE)
    if (is.finite(dev) && dev > 1e-6) {
      warning("Random-effect quadrature not converged (difference ",
              format(signif(dev, 3)), " with a finer rule); increase 'n_points'.",
              call. = FALSE)
    }
  }

  V <- suppressWarnings(stats::vcov(object))
  if (anyNA(V)) {
    if (se_method == "simulation") {
      stop("vcov(object) has missing entries: no draws are possible.", call. = FALSE)
    }
    warning("vcov(object) has missing entries: standard errors are NA.", call. = FALSE)
  }
  est_v <- as.vector(t(est))
  z <- stats::qnorm(1 - (1 - level) / 2)
  jac <- NULL
  if (se_method == "delta") {
    jac <- .brs_num_jacobian(function(th) as.vector(t(quant(th))), theta)
    colnames(jac) <- names(theta)
    se <- sqrt(pmax(rowSums((jac %*% V) * jac), 0))
    se[!is.finite(est_v)] <- NA_real_
    ci <- .brs_marg_wald(est_v, se, z, type, K, contrast = !is.null(pairs))
  } else {
    draws <- .brs_me_rmvnorm(n_sim, theta, V)
    sims <- vapply(seq_len(n_sim), function(r) {
      as.vector(t(quant(stats::setNames(draws[r, ], names(theta)))))
    }, numeric(length(est_v)))
    sims <- matrix(sims, nrow = length(est_v))
    se <- apply(sims, 1L, stats::sd, na.rm = TRUE)
    a <- (1 - level) / 2
    ci <- t(apply(sims, 1L, function(v) {
      if (all(is.na(v))) return(c(NA_real_, NA_real_))
      stats::quantile(v, c(a, 1 - a), na.rm = TRUE, names = FALSE)
    }))
    ci <- list(lower = ci[, 1L], upper = ci[, 2L])
    se[!is.finite(est_v)] <- NA_real_
    ci$lower[!is.finite(est_v)] <- NA_real_
    ci$upper[!is.finite(est_v)] <- NA_real_
  }

  out <- .brs_marg_frame(est, se, ci, type, K, k, pairs, des$newdata)
  attr(out, "type") <- type
  attr(out, "conditional") <- isTRUE(conditional) || !inherits(object, "brsmm")
  attr(out, "level") <- level
  attr(out, "se_method") <- se_method
  attr(out, "n_nodes") <- if (sd_max > 0) length(gh$x) else NA_integer_
  attr(out, "jacobian") <- jac
  class(out) <- c("brs_marginal", "data.frame")
  out
}

# Design matrices of a brs/brsmm fit for `newdata` (NULL: the fit's rows).
# Rows with missing covariates are kept (flag `ok`) so outputs align with newdata.
.brs_marg_design <- function(object, newdata) {
  mixed <- inherits(object, "brsmm")
  mm <- object$model_matrices
  # Variable dispersion: a stored Z (brs fixed dispersion has none)
  has_z <- !is.null(mm$Z)
  if (is.null(newdata)) {
    X <- mm$X
    Z <- if (has_z) mm$Z else NULL
    Xr <- if (mixed) mm$Xr else NULL
  } else {
    if (!is.data.frame(newdata)) {
      stop("'newdata' must be a data.frame.", call. = FALSE)
    }
    if (nrow(newdata) < 1L) stop("'newdata' has no rows.", call. = FALSE)
    build <- function(tt, xlev) {
      tt <- stats::delete.response(tt)
      mf <- stats::model.frame(tt, data = newdata, xlev = xlev,
                               na.action = stats::na.pass)
      M <- stats::model.matrix(tt, mf)
      if (nrow(M) != nrow(newdata)) {
        stop("Could not build the design matrices for 'newdata'.", call. = FALSE)
      }
      M
    }
    X <- build(object$terms$mean, object$xlevels$mean)
    Z <- if (has_z) build(object$terms$precision, object$xlevels$precision) else NULL
    Xr <- if (mixed) build(object$random$re_terms, object$xlevels$random) else NULL
  }
  ok <- stats::complete.cases(X)
  if (!is.null(Z)) ok <- ok & stats::complete.cases(Z)
  if (!is.null(Xr)) ok <- ok & stats::complete.cases(Xr)
  list(
    X = X[ok, , drop = FALSE],
    Z = if (is.null(Z)) NULL else Z[ok, , drop = FALSE],
    Xr = if (is.null(Xr)) NULL else Xr[ok, , drop = FALSE],
    ok = ok, n = length(ok), newdata = newdata
  )
}

# Row pairs (i, j) of a contrast specification as a two-column integer matrix.
.brs_marg_pairs <- function(contrast, n) {
  if (is.null(contrast)) {
    return(NULL)
  }
  pr <- if (is.list(contrast) && !is.data.frame(contrast)) {
    if (!all(lengths(contrast) == 2L)) {
      stop("Each element of 'contrast' must be a pair c(i, j).", call. = FALSE)
    }
    do.call(rbind, lapply(contrast, as.numeric))
  } else if (is.matrix(contrast) || is.data.frame(contrast)) {
    as.matrix(contrast)
  } else if (is.numeric(contrast) && length(contrast) == 2L) {
    matrix(contrast, nrow = 1L)
  } else {
    stop("'contrast' must be a two-column matrix, a list of pairs or a pair c(i, j).",
         call. = FALSE)
  }
  if (ncol(pr) != 2L || !is.numeric(pr) || anyNA(pr) || any(pr != round(pr)) ||
      any(pr < 1) || any(pr > n)) {
    stop("'contrast' must give pairs of row numbers in 1..", n, ".", call. = FALSE)
  }
  storage.mode(pr) <- "integer"
  pr
}

# Largest SD of u = x_r' b over the complete rows (0 without random effects).
.brs_marg_sd_max <- function(object, des, theta) {
  if (is.null(des$Xr) || !nrow(des$Xr)) {
    return(0)
  }
  p <- ncol(des$X)
  q <- if (is.null(des$Z)) 1L else ncol(des$Z)
  q_re <- ncol(des$Xr)
  L <- .brsmm_unpack_chol(unname(theta)[p + q + seq_len(q_re * (q_re + 1L) / 2L)], q_re)
  max(sqrt(rowSums((des$Xr %*% L)^2)))
}

# Nodes z and weights w for E f(Z), Z ~ N(0, 1): composite Gauss-Legendre with
# n_points nodes per panel on [-8.5, 8.5] (tail mass < 2e-17), panels of width
# <= 1 in z and <= 0.5 in u = sd z. `refine` halves the panel width.
.brs_marg_rule <- function(sd_max, n_points, refine = FALSE) {
  if (!is.finite(sd_max) || sd_max <= 0) {
    return(list(x = 0, w = 1))
  }
  zmax <- 8.5
  np <- max(ceiling(2 * zmax), ceiling(2 * zmax * sd_max / 0.5))
  if (refine) np <- 2L * np
  m <- as.integer(n_points)
  # Gauss-Legendre on [-1, 1] (Golub-Welsch)
  kk <- seq_len(m - 1L)
  J <- matrix(0, m, m)
  off <- kk / sqrt(4 * kk^2 - 1)
  J[cbind(kk, kk + 1L)] <- off
  J[cbind(kk + 1L, kk)] <- off
  e <- eigen(J, symmetric = TRUE)
  t <- e$values
  wt <- 2 * e$vectors[1L, ]^2
  br <- seq(-zmax, zmax, length.out = np + 1L)
  half <- diff(br) / 2
  mid <- br[-1L] - half
  x <- as.vector(outer(t, seq_len(np), function(tt, j) mid[j] + tt * half[j]))
  w <- as.vector(outer(wt, seq_len(np), function(ww, j) ww * half[j])) * stats::dnorm(x)
  list(x = x, w = w / sum(w))
}

# Score probabilities of the complete rows at parameters `theta`: quadrature
# average over u ~ N(0, x_r' D x_r) of the cell probabilities at eta + u.
.brs_marg_pmf <- function(theta, des, object, gh) {
  K <- as.integer(object$ncuts)
  n <- nrow(des$X)
  if (n == 0L) {
    return(matrix(numeric(0), 0L, K + 1L))
  }
  p <- ncol(des$X)
  theta <- unname(theta)
  eta_mu <- as.numeric(des$X %*% theta[seq_len(p)])
  if (is.null(des$Z)) {
    eta_phi <- rep(theta[p + 1L], n)
    q <- 1L
  } else {
    q <- ncol(des$Z)
    eta_phi <- as.numeric(des$Z %*% theta[p + seq_len(q)])
  }
  sd_u <- rep(0, n)
  if (!is.null(des$Xr) && length(gh$x) > 1L) {
    q_re <- ncol(des$Xr)
    L <- .brsmm_unpack_chol(theta[p + q + seq_len(q_re * (q_re + 1L) / 2L)], q_re)
    # Var(x_r' b) = |L' x_r|^2
    sd_u <- sqrt(rowSums((des$Xr %*% L)^2))
  }
  ends <- .brs_marg_ends(object, theta)
  m <- length(gh$x)
  eta <- eta_mu + outer(sd_u, gh$x)
  sh <- .brs_shapes_cpp(as.vector(eta), rep(eta_phi, times = m),
                        object$link, object$link_phi, object$repar)
  Pk <- .brs_cell_probs(sh$a, sh$b, ends$left, ends$right, K)
  # Weighted sum over nodes: rows of Pk are (row i, node j), i varying fastest
  P <- vapply(seq_len(K + 1L), function(s) {
    drop(matrix(Pk[, s], n, m) %*% gh$w)
  }, numeric(n))
  P <- matrix(P, nrow = n)
  colnames(P) <- paste0("score_", 0:K)
  P
}

# Cell ends of scores 0..K at parameters `theta`, with the border clamp of the
# likelihood (estimated cutpoints for boundary fits).
.brs_marg_ends <- function(object, theta) {
  K <- as.integer(object$ncuts)
  interval <- .brs_interval_of(object)
  s <- 0:K
  cell <- if (identical(.brs_cut_type_of(object), "boundary")) {
    idx <- .brs_cut_index(object)
    cuts <- .brs_cut_values(theta[idx], K, .brs_cut_anchors(K, interval))$value
    .brs_cell_cut(s, K, interval, cuts)
  } else {
    .brs_cell(s, K, object$lim, interval)
  }
  eps <- .brs_cut_eps
  list(left = pmin(pmax(cell$left, eps), 1 - eps),
       right = pmin(pmax(cell$right, eps), 1 - eps))
}

# Matrix of cell probabilities F(u_s) - F(l_s) (length(a) x (K + 1)); score 0
# and K are the one-sided tails. Interior cells use the tail on the far side of
# the mean, avoiding cancellation; no floor, so partitioning cells sum to 1.
.brs_cell_probs <- function(a, b, left, right, K) {
  N <- length(a)
  P <- matrix(0, N, K + 1L)
  P[, 1L] <- stats::pbeta(right[1L], a, b)
  P[, K + 1L] <- stats::pbeta(left[K + 1L], a, b, lower.tail = FALSE)
  mean_ab <- a / (a + b)
  for (s in seq_len(K - 1L)) {
    l <- left[s + 1L]
    u <- right[s + 1L]
    low <- 0.5 * (l + u) <= mean_ab
    low[is.na(low)] <- TRUE
    v <- numeric(N)
    if (any(low)) {
      v[low] <- stats::pbeta(u, a[low], b[low]) - stats::pbeta(l, a[low], b[low])
    }
    if (any(!low)) {
      v[!low] <- stats::pbeta(l, a[!low], b[!low], lower.tail = FALSE) -
        stats::pbeta(u, a[!low], b[!low], lower.tail = FALSE)
    }
    P[, s + 1L] <- pmax(v, 0)
  }
  P
}

# Jacobian of a vector function by Richardson central differences
# (steps h and h/2, h = 1e-4 max(1, |x_j|)).
.brs_num_jacobian <- function(f, x) {
  f0 <- f(x)
  J <- matrix(NA_real_, length(f0), length(x))
  for (j in seq_along(x)) {
    h <- 1e-4 * max(1, abs(x[j]))
    at <- function(d) {
      xx <- x
      xx[j] <- xx[j] + d
      f(xx)
    }
    d1 <- (at(h) - at(-h)) / (2 * h)
    d2 <- (at(h / 2) - at(-h / 2)) / h
    J[, j] <- (4 * d2 - d1) / 3
  }
  J
}

# Wald intervals: logit scale for probabilities and E[S]/K, natural scale for
# contrasts; a value at the edge of its range gives a degenerate interval.
.brs_marg_wald <- function(est, se, z, type, K, contrast) {
  if (contrast) {
    return(list(lower = est - z * se, upper = est + z * se))
  }
  sc <- if (type == "expected_score") K else 1
  p <- est / sc
  # Values within 1e-10 of an end are degenerate (e.g. P(S >= 0) = 1)
  inside <- is.finite(p) & p > 1e-10 & p < 1 - 1e-10
  lo <- hi <- est
  if (any(inside)) {
    pl <- p[inside]
    lse <- se[inside] / (sc * pl * (1 - pl))
    lo[inside] <- sc * stats::plogis(stats::qlogis(pl) - z * lse)
    hi[inside] <- sc * stats::plogis(stats::qlogis(pl) + z * lse)
  }
  lo[is.na(se)] <- NA_real_
  hi[is.na(se)] <- NA_real_
  list(lower = lo, upper = hi)
}

# Long data frame of the marginal quantities (rows by newdata row or contrast).
.brs_marg_frame <- function(est, se, ci, type, K, k, pairs, newdata) {
  nr <- nrow(est)
  nc <- ncol(est)
  idx <- rep(seq_len(nr), each = nc)
  if (is.null(pairs)) {
    lead <- data.frame(row = idx)
    if (!is.null(newdata)) {
      nd <- newdata[idx, , drop = FALSE]
      rownames(nd) <- NULL
      keep <- setdiff(names(nd), c("row", "score", "k", "estimate", "std.error",
                                    "lower", "upper"))
      lead <- cbind(nd[, keep, drop = FALSE], lead)
    }
  } else {
    lead <- data.frame(
      contrast = paste(pairs[idx, 1L], "-", pairs[idx, 2L]),
      row1 = pairs[idx, 1L], row2 = pairs[idx, 2L],
      stringsAsFactors = FALSE
    )
  }
  extra <- switch(type,
    prob = data.frame(score = rep(0:K, times = nr)),
    exceed = data.frame(k = rep(k, times = nr)),
    expected_score = NULL
  )
  out <- if (is.null(extra)) lead else cbind(lead, extra)
  out$estimate <- as.vector(t(est))
  out$std.error <- se
  out$lower <- ci$lower
  out$upper <- ci$upper
  rownames(out) <- NULL
  out
}
