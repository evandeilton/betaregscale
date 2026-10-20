# ============================================================================ #
# Estimated boundary cutpoints (cutpoints = "boundary")
#
# Score s has the cell (l_s, u_s) of the interval map. The boundary model frees
# the cutpoint c_lo between scores 0 and 1 and the cutpoint c_hi between
# scores K - 1 and K; every other cutpoint keeps the value of the map.
#   K >= 3: c_lo = u_1 plogis(g_lo),  c_hi = l_{K-1} + (1 - l_{K-1}) plogis(g_hi)
#   K  = 2: c_lo = plogis(g_lo),      c_hi = c_lo + (1 - c_lo) plogis(g_hi)
# Both keep the cutpoints ordered for any real (g_lo, g_hi).
# ============================================================================ #

# Admissible values of `cutpoints` (first = default).
.brs_cutpoint_types <- c("fixed", "boundary")

# Partial matching as match.arg(), with a clear message.
.brs_match_cutpoints <- function(cutpoints) {
  if (is.null(cutpoints)) {
    return("fixed")
  }
  tryCatch(match.arg(cutpoints, .brs_cutpoint_types), error = function(e) {
    stop("`cutpoints` must be \"fixed\" or \"boundary\".", call. = FALSE)
  })
}

# Cutpoint type of a fit ("fixed" for objects without the field).
.brs_cut_type_of <- function(object) {
  ct <- object$cutpoints
  if (is.null(ct)) "fixed" else ct
}

# Estimated (c_lo, c_hi) on (0, 1) of a boundary fit; NULL under the fixed map.
.brs_cuts_of <- function(object) {
  if (identical(.brs_cut_type_of(object), "boundary")) object$cut$value else NULL
}

# Fixed-map anchors: c_lo, c_hi and the fixed ends u_1 (upper end of score 1)
# and l_{K-1} (lower end of score K - 1); lim = 0.5, so the cells partition (0, 1).
.brs_cut_anchors <- function(K, interval) {
  cell <- .brs_cell(c(0, 1, K - 1, K), K, 0.5, interval)
  list(c_lo = cell$right[1L], u1 = cell$right[2L],
       l_km1 = cell$left[3L], c_hi = cell$left[4L])
}

# (c_lo, c_hi) and the Jacobian J[k, j] = d c_k / d g_j from the logit-scale g.
.brs_cut_values <- function(g, K, anchors) {
  g <- unname(g)
  p <- stats::plogis(g)
  # plogis'(g) = plogis(g) plogis(-g): no cancellation in 1 - p
  dp <- p * stats::plogis(-g)
  if (K >= 3L) {
    lo <- anchors$u1 * p[1L]
    hi <- anchors$l_km1 + (1 - anchors$l_km1) * p[2L]
    J <- diag(c(anchors$u1 * dp[1L], (1 - anchors$l_km1) * dp[2L]), 2L)
  } else {
    # K = 2: score 1 lies between the two free cutpoints (stick-breaking order)
    lo <- p[1L]
    hi <- lo + (1 - lo) * p[2L]
    J <- matrix(c(dp[1L], stats::plogis(-g[2L]) * dp[1L], 0, (1 - lo) * dp[2L]), 2L, 2L)
  }
  list(value = c(c_lo = lo, c_hi = hi), jacobian = J)
}

# Logit-scale g of given cutpoints (inverse of .brs_cut_values()).
.brs_cut_par <- function(cuts, K, anchors) {
  if (K >= 3L) {
    c(stats::qlogis(cuts[1L] / anchors$u1),
      stats::qlogis((cuts[2L] - anchors$l_km1) / (1 - anchors$l_km1)))
  } else {
    c(stats::qlogis(cuts[1L]), stats::qlogis((cuts[2L] - cuts[1L]) / (1 - cuts[1L])))
  }
}

# Check user cutpoints (c_lo, c_hi) on (0, 1) against the order of the map.
.brs_cut_check_values <- function(cuts, K, interval) {
  if (!is.numeric(cuts) || length(cuts) != 2L || any(!is.finite(cuts))) {
    stop("`cutpoints` must be NULL or two finite numbers c(c_lo, c_hi) on (0, 1).",
         call. = FALSE)
  }
  K <- as.integer(K)
  if (K < 2L) stop("Boundary cutpoints need ncuts >= 2.", call. = FALSE)
  a <- .brs_cut_anchors(K, interval)
  lo_max <- if (K >= 3L) a$u1 else 1
  hi_min <- if (K >= 3L) a$l_km1 else cuts[1L]
  if (!(cuts[1L] > 0 && cuts[1L] < lo_max && cuts[2L] > hi_min && cuts[2L] < 1)) {
    stop("`cutpoints` must satisfy 0 < c_lo < ", if (K >= 3L) format(lo_max) else "c_hi",
         " and ", if (K >= 3L) format(hi_min) else "c_lo", " < c_hi < 1 ",
         "(the cell of each score keeps a positive width).", call. = FALSE)
  }
  as.numeric(cuts)
}

# Scores from values on (0, 1) with the boundary cutpoints `cuts` (the map for
# the interior cutpoints); the coarsening of the boundary model.
.brs_score_from_unit_cut <- function(y, K, interval, cuts) {
  s <- .brs_score_from_unit(y, K, interval)
  if (is.null(cuts)) {
    return(s)
  }
  # Score 0 below c_lo, score K at or above c_hi; between them, the map's
  # border scores move to their neighbours
  s[s == 0 & y >= cuts[1L]] <- 1
  s[s == K & y < cuts[2L]] <- K - 1
  s[y < cuts[1L]] <- 0
  s[y >= cuts[2L]] <- K
  s
}

#' Validate score data for estimated boundary cutpoints
#'
#' @param Y Response matrix from \code{.extract_response()}.
#' @param K,lim,interval Grid of the fit.
#' @param require_all Stop when a score next to a free cutpoint is never
#'   observed (the cutpoint is then not identified).
#' @return A list with the scores, row indices of the four border roles,
#'   the fixed-map anchors and starting values \code{start} (the map's g).
#' @keywords internal
#' @noRd
.brs_cut_setup <- function(Y, K, lim, interval, require_all = TRUE) {
  K <- as.integer(K)
  if (K < 2L) {
    stop("cutpoints = \"boundary\" needs ncuts >= 2.", call. = FALSE)
  }
  if (abs(lim - 0.5) > sqrt(.Machine$double.eps)) {
    stop("cutpoints = \"boundary\" needs lim = 0.5 (cells that partition (0, 1)); ",
         "found lim = ", format(lim), ".", call. = FALSE)
  }
  y <- as.numeric(Y[, "y"])
  delta <- as.integer(Y[, "delta"])
  n_exact <- sum(delta == 0L)
  if (n_exact > 0L) {
    stop("cutpoints = \"boundary\" needs score data: ", n_exact, " exact ",
         "observation(s) (delta = 0) found. Use cutpoints = \"fixed\".", call. = FALSE)
  }
  # Every row must be exactly the cell of its score under the map
  is_score <- is.finite(y) & abs(y - round(y)) < 1e-8 & y >= 0 & y <= K
  s <- round(y)
  ok <- is_score
  if (any(is_score)) {
    ref <- .brs_check_core(s[is_score], ncuts = K, lim = 0.5, delta = NULL,
                           interval = interval)
    ok[is_score] <- ref[, "delta"] == delta[is_score] &
      abs(ref[, "left"] - Y[is_score, "left"]) < 1e-9 &
      abs(ref[, "right"] - Y[is_score, "right"]) < 1e-9
  }
  if (!all(ok)) {
    bad <- which(!ok)
    stop("cutpoints = \"boundary\" needs every row to be the cell of a score 0..",
         K, ": ", length(bad), " row(s) carry analyst-supplied censoring or ",
         "thresholds (first: row ", bad[1L], "). Use cutpoints = \"fixed\".",
         call. = FALSE)
  }
  need <- if (K == 2L) 0:2 else c(0L, 1L, K - 1L, K)
  miss <- setdiff(need, s)
  if (require_all && length(miss)) {
    stop("cutpoints = \"boundary\": score(s) ", paste(miss, collapse = ", "),
         " not observed, so the boundary cutpoint next to them is not identified ",
         "(the likelihood keeps increasing towards a bound). Use cutpoints = \"fixed\".",
         call. = FALSE)
  }
  anchors <- .brs_cut_anchors(K, interval)
  list(
    K = K, interval = interval, s = s, delta = delta,
    left = as.numeric(Y[, "left"]), right = as.numeric(Y[, "right"]),
    yt = as.numeric(Y[, "yt"]),
    i0 = which(s == 0), i1 = which(s == 1), ikm1 = which(s == K - 1), ik = which(s == K),
    rows = which(s %in% need),
    anchors = anchors,
    fixed = c(c_lo = anchors$c_lo, c_hi = anchors$c_hi),
    start = .brs_cut_par(c(anchors$c_lo, anchors$c_hi), K, anchors)
  )
}

# Cell endpoints of all rows at cutpoints `cv` (from .brs_cut_values()).
.brs_cut_endpoints <- function(cv, cs) {
  eps <- 1e-5
  # Same border clamp as .brs_check_core()
  cl <- min(max(cv$value[1L], eps), 1 - eps)
  ch <- min(max(cv$value[2L], eps), 1 - eps)
  left <- cs$left
  right <- cs$right
  right[cs$i0] <- cl
  left[cs$i1] <- cl
  right[cs$ikm1] <- ch
  left[cs$ik] <- ch
  list(left = left, right = right)
}

# Beta shapes as the compiled likelihood builds them (src/brs_common.h).
.brs_shapes_cpp <- function(eta_mu, eta_phi, link, link_phi, repar) {
  mu <- .clamp_mu_by_repar(apply_inv_link(eta_mu, link), repar)
  phi <- .clamp_phi_by_repar(apply_inv_link(eta_phi, link_phi), repar)
  if (repar == 0L) {
    a <- mu
    b <- phi
  } else if (repar == 1L) {
    a <- mu * phi
    b <- (1 - mu) * phi
  } else {
    r <- (1 - phi) / phi
    a <- mu * r
    b <- (1 - mu) * r
  }
  list(a = pmin(pmax(a, 1e-12), 1e8), b = pmin(pmax(b, 1e-12), 1e8))
}

# Analytic d loglik / d (c_lo, c_hi): +-f(c) / P_i for the rows whose cell
# ends (upper: +, lower: -) at c. `sh` holds the shapes of rows cs$rows.
.brs_cut_grad_c <- function(cv, ends, cs, sh) {
  eps <- 1e-5
  r <- cs$rows
  ll <- .brs_obs_loglik(cs$delta[r], ends$left[r], ends$right[r], cs$yt[r],
                        sh$a, sh$b)
  term <- function(idx, cval) {
    if (!length(idx) || cval <= eps || cval >= 1 - eps) {
      return(0)  # clamped cutpoint: the likelihood is flat in it
    }
    j <- match(idx, r)
    li <- ll[j]
    keep <- li > -1e6  # penalised rows are constant in c
    v <- exp(stats::dbeta(cval, sh$a[j], sh$b[j], log = TRUE) - li)
    sum(v[keep & is.finite(v)])
  }
  c_lo <- cv$value[1L]
  c_hi <- cv$value[2L]
  c(term(cs$i0, c_lo) - term(cs$i1, c_lo),
    term(cs$ikm1, c_hi) - term(cs$ik, c_hi))
}

# Richardson central difference of f at x along coordinate j (steps h, h/2).
.brs_rich_diff <- function(f, x, j, h = 1e-3 * max(1, abs(x[j]))) {
  at <- function(d) {
    xx <- x
    xx[j] <- xx[j] + d
    f(xx)
  }
  d1 <- (at(h) - at(-h)) / (2 * h)
  d2 <- (at(h / 2) - at(-h / 2)) / h
  (4 * d2 - d1) / 3
}

#' Objective, gradient and Hessian of a model with boundary cutpoints
#'
#' @param cs Setup from \code{.brs_cut_setup()}.
#' @param k0 Number of parameters before the two cutpoint parameters.
#' @param ll0,gr0,hess0 Log-likelihood, gradient and Hessian in the first
#'   \code{k0} parameters given the endpoints: \code{f(par0, left, right)}.
#' @param shapes Function \code{(par0, rows)} returning the beta shapes of
#'   \code{rows}, or \code{NULL} to differentiate in g numerically (mixed models).
#' @return A list of functions \code{fn}, \code{gr}, \code{hess},
#'   \code{ends} (endpoints at a full parameter vector).
#' @keywords internal
#' @noRd
.brs_cut_problem <- function(cs, k0, ll0, gr0, hess0, shapes = NULL) {
  idx0 <- seq_len(k0)
  idx_g <- k0 + 1:2
  at <- function(g) {
    cv <- .brs_cut_values(g, cs$K, cs$anchors)
    list(cv = cv, ends = .brs_cut_endpoints(cv, cs))
  }
  fn <- function(par) {
    e <- at(par[idx_g])$ends
    ll0(par[idx0], e$left, e$right)
  }
  gr <- function(par) {
    par0 <- par[idx0]
    a <- at(par[idx_g])
    g0 <- gr0(par0, a$ends$left, a$ends$right)
    gg <- if (is.null(shapes)) {
      # Mixed models: Richardson differences of the marginal likelihood in g
      vapply(1:2, function(j) .brs_rich_diff(fn, par, k0 + j), numeric(1))
    } else {
      sh <- shapes(par0, cs$rows)
      as.numeric(crossprod(a$cv$jacobian, .brs_cut_grad_c(a$cv, a$ends, cs, sh)))
    }
    c(g0, gg)
  }
  hess <- function(par) {
    par0 <- par[idx0]
    g <- par[idx_g]
    e <- at(g)$ends
    H <- matrix(0, k0 + 2L, k0 + 2L)
    H[idx0, idx0] <- hess0(par0, e$left, e$right)
    # Cross block: derivative of the first-block gradient in g
    Jx <- numDeriv::jacobian(function(gg) {
      ee <- at(gg)$ends
      gr0(par0, ee$left, ee$right)
    }, g)
    H[idx0, idx_g] <- Jx
    H[idx_g, idx0] <- t(Jx)
    # Cutpoint block: Jacobian of the analytic gradient, or Hessian of the
    # marginal likelihood in g (mixed models)
    Hg <- if (is.null(shapes)) {
      numDeriv::hessian(function(gg) fn(c(par0, gg)), g)
    } else {
      numDeriv::jacobian(function(gg) gr(c(par0, gg))[idx_g], g)
    }
    H[idx_g, idx_g] <- 0.5 * (Hg + t(Hg))
    H
  }
  list(fn = fn, gr = gr, hess = hess, ends = function(par) at(par[idx_g])$ends)
}

# Boundary-cutpoint problem of a brs fit: compiled likelihood, gradient and
# Hessian in (beta, gamma) at the moved cell ends, analytic gradient in g.
# Z = NULL: fixed dispersion (one coefficient after beta).
.brs_cut_problem_brs <- function(cs, X, Z, Y, delta, link, link_phi, repar) {
  lc_mu <- link_to_code(link)
  lc_phi <- link_to_code(link_phi)
  yt <- as.numeric(Y[, "yt"])
  p <- ncol(X)
  if (is.null(Z)) {
    .brs_cut_problem(
      cs, k0 = p + 1L,
      ll0 = function(par0, l, r) {
        .brs_loglik_fixed_cpp(par0, X, l, r, yt, delta, lc_mu, lc_phi, repar)
      },
      gr0 = function(par0, l, r) {
        .brs_grad_fixed_cpp(par0, X, l, r, yt, delta, lc_mu, lc_phi, repar)
      },
      hess0 = function(par0, l, r) {
        .brs_hessian_fixed_cpp(par0, X, l, r, yt, delta, lc_mu, lc_phi, repar)
      },
      shapes = function(par0, rows) {
        .brs_shapes_cpp(X[rows, , drop = FALSE] %*% par0[seq_len(p)], par0[p + 1L],
                        link, link_phi, repar)
      }
    )
  } else {
    q <- ncol(Z)
    .brs_cut_problem(
      cs, k0 = p + q,
      ll0 = function(par0, l, r) {
        .brs_loglik_variable_cpp(par0, X, Z, l, r, yt, delta, lc_mu, lc_phi, repar)
      },
      gr0 = function(par0, l, r) {
        .brs_grad_variable_cpp(par0, X, Z, l, r, yt, delta, lc_mu, lc_phi, repar)
      },
      hess0 = function(par0, l, r) {
        .brs_hessian_variable_cpp(par0, X, Z, l, r, yt, delta, lc_mu, lc_phi, repar)
      },
      shapes = function(par0, rows) {
        .brs_shapes_cpp(X[rows, , drop = FALSE] %*% par0[seq_len(p)],
                        Z[rows, , drop = FALSE] %*% par0[p + seq_len(q)],
                        link, link_phi, repar)
      }
    )
  }
}

# Response matrix of a boundary fit with the fixed-map cells restored.
.brs_cut_fixed_Y <- function(object) {
  Y <- object$Y
  ref <- .brs_check_core(round(Y[, "y"]), ncuts = object$ncuts, lim = 0.5,
                         delta = NULL, interval = .brs_interval_of(object))
  Y[, "left"] <- ref[, "left"]
  Y[, "right"] <- ref[, "right"]
  Y
}

# Log-likelihood and gradient of a brs fit with boundary cutpoints as
# functions of the full parameter vector (checks and profiling).
.brs_cut_objective <- function(object) {
  Y <- .brs_cut_fixed_Y(object)
  cs <- .brs_cut_setup(Y, object$ncuts, object$lim, .brs_interval_of(object),
                       require_all = FALSE)
  .brs_cut_problem_brs(cs, object$model_matrices$X, object$model_matrices$Z, Y,
                       as.integer(Y[, "delta"]), object$link, object$link_phi,
                       object$repar)
}

# The `cut` component of a boundary fit: estimates on the logit and unit scales.
# A logit-scale estimate beyond +-10 means a border cell is collapsing (advisory).
.brs_cut_component <- function(cs, g_hat) {
  far <- abs(g_hat) > 10 | !is.finite(g_hat)
  if (any(far)) {
    .brs_advisory(paste0(
      "Boundary cutpoint(s) ", paste(c("c_lo", "c_hi")[far], collapse = ", "),
      " at the edge of their range (logit-scale estimate ",
      paste(format(g_hat[far], digits = 3), collapse = ", "), "): a cell next to ",
      "it is nearly empty and the cutpoint is not identified; use ",
      "cutpoints = \"fixed\" or merge sparse categories."
    ))
  }
  cv <- .brs_cut_values(g_hat, cs$K, cs$anchors)
  list(
    par = stats::setNames(unname(as.numeric(g_hat)), c("(cut_lo)", "(cut_hi)")),
    value = cv$value,
    jacobian = cv$jacobian,
    fixed = cs$fixed,
    start = cs$start,
    K = cs$K,
    interval = cs$interval
  )
}

# Table of the boundary cutpoints: logit-scale estimate, the cutpoint on (0, 1)
# and on the latent score scale, with delta-method standard errors.
.brs_cut_table <- function(object, V) {
  cut <- object$cut
  if (is.null(cut)) {
    return(NULL)
  }
  idx <- match(c("(cut_lo)", "(cut_hi)"), names(object$par))
  Vg <- V[idx, idx, drop = FALSE]
  J <- cut$jacobian
  Vc <- J %*% Vg %*% t(J)
  K <- cut$K
  iv <- if (is.null(cut$interval)) .brs_interval_of(object) else cut$interval
  # Latent score scale of the map: K c (mid), (K + 1) c (right), (K + 1) c - 1 (left)
  sc <- if (identical(iv, "mid")) K else K + 1
  tab <- cbind(
    Estimate = cut$par,
    `Std. Error` = sqrt(diag(Vg)),
    cutpoint = cut$value,
    `SE(cut)` = sqrt(diag(Vc)),
    score = .brs_latent_score(cut$value, K, iv),
    `SE(score)` = sc * sqrt(diag(Vc)),
    fixed = .brs_latent_score(cut$fixed, K, iv)
  )
  rownames(tab) <- c(paste0("c_lo (0|1)"), paste0("c_hi (", K - 1L, "|", K, ")"))
  tab
}

# Print the cutpoint table of a summary.
.brs_print_cut_table <- function(tab, digits) {
  if (is.null(tab)) {
    return(invisible(NULL))
  }
  cat("Boundary cutpoints (logit scale; on (0, 1) and on the score scale, ",
      "with delta-method SEs; 'fixed' = the map's value):\n", sep = "")
  print(round(tab, digits))
  cat("\n")
  invisible(NULL)
}

# Endpoints of a response matrix under the estimated cutpoints of a fit
# (held-out data in brs_cv()); unchanged under the fixed map.
.brs_apply_fit_cuts <- function(Y, object) {
  cuts <- .brs_cuts_of(object)
  if (is.null(cuts)) {
    return(Y)
  }
  cs <- .brs_cut_setup(Y, object$ncuts, object$lim, .brs_interval_of(object),
                       require_all = FALSE)
  e <- .brs_cut_endpoints(list(value = cuts), cs)
  Y[, "left"] <- e$left
  Y[, "right"] <- e$right
  Y
}

# Positions of the cutpoint parameters in coef(); stops for fixed-map fits.
.brs_cut_index <- function(object) {
  idx <- match(c("(cut_lo)", "(cut_hi)"), names(object$par))
  if (anyNA(idx)) {
    stop("The model has no estimated cutpoints (cutpoints = \"fixed\").", call. = FALSE)
  }
  idx
}

# One print() line with the estimated boundary cutpoints.
.brs_print_cuts <- function(x, digits) {
  cuts <- .brs_cuts_of(x)
  if (is.null(cuts)) {
    return(invisible(NULL))
  }
  K <- x$cut$K
  sc <- .brs_latent_score(cuts, K, .brs_interval_of(x))
  cat("Boundary cutpoints (0|1, ", K - 1L, "|", K, "): ",
      paste(formatC(cuts, digits = digits, format = "f"), collapse = ", "),
      " on (0, 1); ", paste(formatC(sc, digits = digits, format = "f"), collapse = ", "),
      " on the score scale\n\n", sep = "")
  invisible(NULL)
}
